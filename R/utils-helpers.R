# ===================================================================
# INTERNAL HELPER FUNCTIONS
# Support for `find_cutpoint()` and `find_cutpoint_number()`.
# ===================================================================

#' Objective Function for Genetic Algorithm
#'
#' @description
#' Called by `rgenoud::genoud`. Calculates fitness (log-likelihood,
#' log-rank stat, etc.) for a given set of cut-points.
#' `rgenoud` always maximises, so all criteria are framed as such.
#'
#' @param params Numeric vector. First `numcut` are cut-points;
#'   remaining (if `criterion = "loglik"`) are initial betas.
#' @param time,censor,target,confound Survival data and predictor.
#' @param numcut Number of cut-points.
#' @param gap Minimum distance between cut-points.
#' @param nmin Minimum observations per group.
#' @param criterion `"logrank"`, `"p_value"`, `"hazard_ratio"` or `"loglik"`.
#' @param loglik0 Log-likelihood of the null model (for `"p_value"`).
#'
#' @return Single numeric fitness value.
#'
#' @references
#' Cox, D. R. (1972). Regression Models and Life-Tables.
#' *Journal of the Royal Statistical Society: Series B*, **34**(2),
#' 187–202. \doi{10.1111/j.2517-6161.1972.tb00899.x}
#'
#' Mantel, N. (1966). Evaluation of survival data and two new rank order
#' statistics arising in its consideration. *Cancer Chemotherapy Reports*,
#' **50**(3). <https://pubmed.ncbi.nlm.nih.gov/5910392/>
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.3}  `na.rm = TRUE` in `quantile`.
#' @srrstats {G2.4a} `make.names()` sanitises covariate names.
#' @srrstats {G2.13} Invalid gaps return `-Inf`.
#' @srrstats {G5.0}  Edge cases (zero gap, infinite domain) return `-Inf`.
#' @srrstats {G5.1}  Failure modes documented.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.obj <- function(params, time, censor, target, confound, numcut, gap, nmin,
                 criterion, loglik0 = NA_real_) {
  cutoff <- params[1:numcut]
  if (numcut > 1 && min(diff(sort(cutoff))) < gap) {
    return(-Inf)
  }
  cut_design <- factor(cut(target,
                           breaks = c(-Inf, sort(cutoff), Inf),
                           labels = FALSE, include.lowest = TRUE
  ))
  if (length(table(cut_design)) != (numcut + 1)) {
    return(-Inf)
  }
  if (min(table(cut_design)) < nmin) {
    return(-Inf)
  }

  data_for_fit <- data.frame(
    time = time, censor = censor,
    cut_design = cut_design
  )
  if (!is.null(confound) && ncol(confound) > 0) {
    colnames(confound) <- make.names(colnames(confound))
    data_for_fit <- cbind(data_for_fit, confound)
  }
  formula_str <- "survival::Surv(time, censor) ~ cut_design"
  if (!is.null(confound) && ncol(confound) > 0) {
    formula_str <- paste(formula_str, "+",
                         paste(colnames(confound), collapse = " + ")
    )
  }
  fit_formula <- as.formula(formula_str)

  stat_value <- switch(
    criterion,
    "logrank" = {
      fit <- tryCatch(survival::survdiff(fit_formula, data = data_for_fit),
                      error = function(e) NULL
      )
      if (is.null(fit)) -Inf else fit$chisq
    },
    "hazard_ratio" = {
      fit <- tryCatch(survival::coxph(fit_formula, data = data_for_fit),
                      error = function(e) NULL
      )
      if (is.null(fit) || is.null(fit$coefficients) ||
          !("cut_design2" %in% names(fit$coefficients))) {
        -Inf
      } else {
        exp(fit$coefficients["cut_design2"])
      }
    },
    "p_value" = {
      fit <- tryCatch(survival::coxph(fit_formula, data = data_for_fit),
                      error = function(e) NULL
      )
      if (is.null(fit) || is.null(fit$loglik)) {
        return(-Inf)
      }
      if (!is.na(loglik0)) {
        loglik1 <- fit$loglik[2]
        statistic <- -2 * (loglik0 - loglik1)
        df <- numcut
        pval <- stats::pchisq(statistic, df, lower.tail = FALSE)
        return(1 - pval)
      } else {
        sfit <- tryCatch(summary(fit), error = function(e) NULL)
        if (is.null(sfit) || is.null(sfit$logtest)) {
          -Inf
        } else {
          1 - sfit$logtest["pvalue"]
        }
      }
    },
    "loglik" = {
      beta_init <- params[(numcut + 1):length(params)]
      fit <- tryCatch(
        survival::coxph(fit_formula,
                        data = data_for_fit,
                        init = beta_init, iter.max = 0
        ),
        error = function(e) NULL
      )
      if (is.null(fit) || is.null(fit$loglik)) -Inf else fit$loglik[2]
    },
    -Inf
  )
  return(stat_value)
}

#' Wrapper for the rgenoud Genetic Algorithm
#'
#' @description
#' Sets up and runs `rgenoud::genoud` to find optimal cut-points.
#'
#' @param target Continuous predictor vector.
#' @param numcut Number of cut-points.
#' @param criterion Optimisation criterion.
#' @param ... Additional arguments passed to `.obj` and `genoud`.
#'
#' @return A `genoud` object (or `NULL` on failure).
#'
#' @references
#' Mebane Jr, W. R., & Sekhon, J. S. (2011). Genetic Optimization Using
#' Derivatives: The rgenoud Package for R. *Journal of Statistical Software*,
#' **42**, 1–26. \doi{10.18637/jss.v042.i11}
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.13} `requireNamespace("rgenoud")` check.
#' @srrstats {G5.0}  Edge cases (zero gap, infinite range) return `NULL`.
#' @srrstats {G5.1}  Failure modes documented.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.run_genetic_search <- function(target, numcut, time, censor, confound, nmin,
                                criterion, numgen = 15, gap = NULL,
                                print.level = 0, ...) {
  num_confound_vars <- if (is.null(confound) || all(is.na(confound))) {
    0
  } else {
    ncol(confound)
  }
  optimizing_betas <- (criterion == "loglik")
  nvars <- if (optimizing_betas) {
    numcut + numcut + num_confound_vars
  } else {
    numcut
  }

  domain_cuts <- matrix(rep(range(target, na.rm = TRUE), numcut),
                        ncol = 2, byrow = TRUE
  )
  if (any(is.infinite(domain_cuts))) {
    return(NULL)
  }

  domain <- if (optimizing_betas) {
    domain_betas <- matrix(rep(c(-5, 5), nvars - numcut),
                           ncol = 2,
                           byrow = TRUE
    )
    rbind(domain_cuts, domain_betas)
  } else {
    domain_cuts
  }

  initial_cuts <- stats::quantile(target,
                                  probs = seq(0, 1, length.out = numcut + 2),
                                  na.rm = TRUE
  )[2:(numcut + 1)]
  initial_values <- if (optimizing_betas) {
    initial_betas <- rep(0, nvars - numcut)
    c(initial_cuts, initial_betas)
  } else {
    initial_cuts
  }
  if (any(is.na(initial_values))) {
    return(NULL)
  }

  if (is.null(gap)) {
    gap <- stats::quantile(sort(diff(sort(unique(na.omit(target))))),
                           probs = 0.5, na.rm = TRUE
    )
    if (is.na(gap) || gap == 0) gap <- 1e-4
  }

  loglik0 <- NA_real_
  if (criterion == "p_value") {
    null_formula_str <- "survival::Surv(time, censor) ~ 1"
    data_for_null_fit <- data.frame(time = time, censor = censor)
    if (!is.null(confound) && ncol(confound) > 0) {
      confound_null <- confound
      colnames(confound_null) <- make.names(colnames(confound_null))
      data_for_null_fit <- cbind(data_for_null_fit, confound_null)
      null_formula_str <- paste(null_formula_str, "+",
                                paste(colnames(confound_null),
                                      collapse = " + "
                                )
      )
    }
    null_fit <- tryCatch(
      survival::coxph(as.formula(null_formula_str), data = data_for_null_fit),
      error = function(e) NULL
    )
    if (!is.null(null_fit) && !is.null(null_fit$loglik)) {
      loglik0 <- null_fit$loglik[length(null_fit$loglik)]
    }
  }

  if (!requireNamespace("rgenoud", quietly = TRUE)) {
    stop("The 'rgenoud' package is required. Please install it.", call. = FALSE)
  }

  optim_result <- tryCatch(rgenoud::genoud(
    fn = .obj,
    nvars = nvars,
    max = TRUE,
    pop.size = 100,
    max.generations = numgen,
    wait.generations = 10,
    hard.generation.limit = TRUE,
    starting.values = initial_values,
    Domains = domain,
    print.level = print.level,
    time = time,
    censor = censor,
    target = target,
    confound = confound,
    numcut = numcut,
    gap = gap,
    nmin = nmin,
    criterion = criterion,
    loglik0 = loglik0
  ), error = function(e) NULL)

  return(optim_result)
}

#' Calculate Information Criterion (AIC, AICc, or BIC)
#'
#' @description
#' Computes AIC, AICc, or BIC from a `coxph` object or log-likelihood.
#'
#' @param model Fitted `coxph` object or list with `loglik` component.
#' @param k Number of parameters (degrees of freedom).
#' @param n Sample size.
#' @param criterion `"AIC"`, `"AICc"` or `"BIC"`.
#'
#' @return Single numeric IC value (or `NA` on failure).
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' *IEEE Transactions on Automatic Control*, **19**(6), 716–723.
#' \doi{10.1109/TAC.1974.1100705}
#'
#' Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model
#' selection in small samples. *Biometrika*, **76**(2), 297–307.
#' \doi{10.1093/biomet/76.2.297}
#'
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' *The Annals of Statistics*, **6**(2), 461–464.
#' \doi{10.1214/aos/1176344136}
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.1} Implements AIC, AICc, BIC.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.calc_ic <- function(model, k, n, criterion) {
  if (is.null(model) || !is.list(model)) {
    return(NA_real_)
  }
  logL <- model$loglik[2]
  if (is.null(logL) || is.na(logL)) {
    return(NA_real_)
  }
  if (criterion == "BIC") {
    return(-2 * logL + k * log(n))
  } else if (criterion == "AICc") {
    if ((n - k - 1) <= 0) {
      return(NA_real_)
    }
    aic <- -2 * logL + 2 * k
    aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
    return(aicc)
  } else { # "AIC"
    return(-2 * logL + 2 * k)
  }
}

# ===================================================================
# Validation & Prep Helpers for find_cutpoint()
# ===================================================================

#' Validate inputs for find_cutpoint()
#'
#' @description
#' Centralised validation for `find_cutpoint()` arguments.
#'
#' @param method,criterion,num_cuts,covariates,predictor,outcome_time,
#'   outcome_event,data Arguments from `find_cutpoint()`.
#'
#' @return `invisible(TRUE)` on success; aborts on failure.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.0}  Validates `method`, `criterion`, `num_cuts`.
#' @srrstats {G2.1}  Type checks.
#' @srrstats {G2.4b} `match.arg()` for controlled vocabularies.
#' @srrstats {G2.4c} `requireNamespace()` for optional `rgenoud`.
#' @srrstats {G2.9}  Column-name existence.
#' @srrstats {G2.13} `cli_abort()` for invalid inputs.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.validate_find_cutpoint_inputs <- function(data, predictor, outcome_time,
                                           outcome_event, num_cuts, method,
                                           criterion, covariates) {
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion,
                         choices = c("logrank", "hazard_ratio", "p_value")
  )
  if (!is.numeric(num_cuts) || num_cuts < 0 || num_cuts != round(num_cuts)) {
    stop("num_cuts must be a non-negative integer.", call. = FALSE)
  }
  if (criterion == "hazard_ratio" && num_cuts > 1) {
    stop("'hazard_ratio' is only supported for num_cuts = 1.",
         call. = FALSE
    )
  }

  if (method == "genetic" && !requireNamespace("rgenoud", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "'genetic' method requires 'rgenoud'.",
        "i" = "Install with: install.packages(\"rgenoud\")",
        "i" = "Or use `method = \"systematic\"` (for num_cuts <= 2)."
      )
    )
  }
  if (method == "systematic" && !num_cuts %in% c(1, 2)) {
    stop("Systematic search only supports num_cuts = 1 or 2.",
         call. = FALSE
    )
  }

  if (is.null(predictor)) {
    stop("A 'predictor' variable must be specified.", call. = FALSE)
  }
  if (is.null(outcome_time) || is.null(outcome_event)) {
    stop("Both 'outcome_time' and 'outcome_event' are required.",
         call. = FALSE
    )
  }

  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(
      paste0(
        "Missing required columns: '",
        paste(missing_vars, collapse = "', '"), "'"
      ),
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Prepare data for cut-point analysis
#'
#' @description
#' Subsets, removes incomplete cases, and renames columns.
#'
#' @param data,predictor,outcome_time,outcome_event,covariates Arguments
#'   from `find_cutpoint()`.
#'
#' @return Cleaned `userdata` data.frame.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.3}  `na.omit()` with explicit `NA` handling.
#' @srrstats {G2.3a} Removes incomplete cases.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.prepare_cutpoint_data <- function(data, predictor, outcome_time,
                                   outcome_event, covariates) {
  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  userdata <- data[, unique(required_vars), drop = FALSE]
  userdata <- stats::na.omit(userdata)
  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"
  return(userdata)
}

#' Validate data conditions for cut-point analysis
#'
#' @description
#' Checks post-cleaning: non-negative time, valid 0/1 event,
#' sufficient data, non-constant predictor.
#'
#' @param userdata Cleaned data from `.prepare_cutpoint_data()`.
#' @param nmin,num_cuts,quiet,outcome_event Args from `find_cutpoint()`.
#'
#' @return List with `valid` (logical) and `nmin_abs` (int).
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.1}  Event column must be numeric 0/1.
#' @srrstats {G2.13} `cli_abort()` for invalid event data.
#' @srrstats {G5.0}  Edge cases (constant predictor, insufficient data).
#' @srrstats {G5.4a} Checks `n < nmin * (num_cuts + 1)`.
#' @srrstats {G5.4b} Checks `unique(factor) <= num_cuts`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.validate_data_conditions <- function(userdata, nmin, num_cuts,
                                      outcome_event, quiet) {
  if (any(userdata$time < 0, na.rm = TRUE)) {
    stop("Time variable must be non-negative.", call. = FALSE)
  }

  event_col <- userdata$event
  if (!is.numeric(event_col)) {
    cli::cli_abort(c(
      "Event column ({.arg {outcome_event}}) must be numeric.",
      "i" = "Detected: {.cls {class(event_col)}}.",
      "i" = "Please convert to 0 (censored) and 1 (event)."
    ))
  }
  valid_events <- unique(event_col)
  if (!all(valid_events %in% c(0, 1))) {
    invalid_vals <- sort(valid_events[!(valid_events %in% c(0, 1))])
    cli::cli_abort(c(
      "Event column ({.arg {outcome_event}}) must contain only 0 and 1.",
      "i" = "Found invalid value(s): {.val {invalid_vals}}"
    ))
  }

  if (nmin > 0 && nmin < 1) {
    nmin_abs <- floor(nmin * nrow(userdata))
    if (!quiet) {
      cli::cli_alert_info(paste(
        "nmin {nmin} is a proportion.",
        "Min. group size set to {nmin_abs}."
      ))
    }
  } else if (nmin >= 1) {
    nmin_abs <- floor(nmin)
  } else {
    stop("'nmin' must be a non-negative number.", call. = FALSE)
  }

  if (nrow(userdata) < nmin_abs * (num_cuts + 1)) {
    if (!quiet) {
      cli::cli_inform(paste(
        "Not enough data ({nrow(userdata)}) for nmin",
        "({nmin_abs}) and {num_cuts} cut(s). Returning NA."
      ))
    }
    return(list(valid = FALSE, nmin_abs = nmin_abs))
  }

  if (length(unique(userdata$factor)) <= num_cuts) {
    if (!quiet) {
      cli::cli_inform(paste(
        "Predictor has too few unique values",
        "({length(unique(userdata$factor))}) for {num_cuts} cut(s).",
        "Returning NA."
      ))
    }
    return(list(valid = FALSE, nmin_abs = nmin_abs))
  }

  return(list(valid = TRUE, nmin_abs = nmin_abs))
}

#' Infix operator for NULL default
#'
#' @param a The value to check.
#' @param b The default value if `a` is `NULL`.
#'
#' @return `a` if not `NULL`, otherwise `b`.
#'
#' @examples
#' x <- NULL
#' y <- 5
#' x %||% y # Returns 5
#'
#' z <- 10
#' z %||% y # Returns 10
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.0} Input validation.
#' @srrstats {G2.1} Type checking.
#'
#' @name or-operator
#' @export
`%||%` <- function(a, b) if (is.null(a)) b else a
