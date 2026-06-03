# ===================================================================
# INTERNAL ENGINE: GENETIC SEARCH (rgenoud)
# Handles evolutionary threshold discovery for multi-dimensional cuts
# ===================================================================

#' Internal helper: Objective Function for Genetic Algorithm
#'
#' @description
#' Called by `rgenoud::genoud`. Calculates fitness (log-likelihood,
#' log-rank stat, etc.) for a given set of cut-points.
#' `rgenoud` always maximises, so all criteria are framed as such.
#'
#' @param params Numeric vector. First `numcut` are cut-points; remaining are betas.
#' @param time Survival time vector.
#' @param censor Survival event vector.
#' @param target Continuous predictor vector.
#' @param confound Optional covariate data frame.
#' @param numcut Number of cut-points.
#' @param gap Minimum distance between cut-points.
#' @param nmin Minimum observations per group.
#' @param criterion "logrank", "p_value", "hazard_ratio" or "loglik".
#' @param loglik0 Log-likelihood of the null model.
#' @param cache Optional environment to cache and retrieve evaluations.
#' @param base_df Pre-allocated data.frame template for speed.
#' @param precompiled_formula Pre-parsed survival formula.
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
#' @srrstats {G2.3} `na.rm = TRUE` used in quantile calculations.
#' @srrstats {G2.4a} `make.names()` sanitises covariate names.
#' @srrstats {G2.4d} Explicit conversion to factor used for group creation.
#' @srrstats {G2.13} Invalid gaps or insufficient N return `-Inf`.
#' @srrstats {G5.8} Edge cases (zero gap, infinite domain) return `-Inf`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival survdiff coxph
#' @importFrom stats pchisq
#' @noRd
.obj <- function(params, time, censor, target, confound, numcut, gap, nmin,
                 criterion, loglik0 = NA_real_, cache = NULL,
                 base_df = NULL, precompiled_formula = NULL) {
  if (!is.null(cache)) {
    key <- paste(round(params, 6), collapse = "_")
    if (exists(key, envir = cache, inherits = FALSE)) return(cache[[key]])
  }

  if (length(unique(time)) <= 1) return(-Inf)
  cutoff <- params[1:numcut]
  sorted_cuts <- sort(cutoff)

  if (numcut > 1 && min(diff(sorted_cuts)) < gap) return(-Inf)
  if (sorted_cuts[1] <= min(target, na.rm = TRUE) || sorted_cuts[numcut] >= max(target, na.rm = TRUE)) return(-Inf)

  cut_design <- as.factor(findInterval(target, sorted_cuts, left.open = TRUE) + 1L)
  if (length(table(cut_design)) != (numcut + 1) || min(table(cut_design)) < nmin) return(-Inf)

  data_for_fit <- base_df
  data_for_fit$cut_design <- cut_design

  stat_value <- switch(criterion,
                       "logrank" = {
                         if (is.null(confound) || ncol(confound) == 0) {
                           fit <- tryCatch(survival::survdiff(precompiled_formula, data = data_for_fit), error = function(e) NULL)
                           if (is.null(fit)) -Inf else fit$chisq
                         } else {
                           fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
                           if (is.null(fit) || is.null(fit$score)) -Inf else fit$score
                         }
                       },
                       "hazard_ratio" = {
                         fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
                         if (is.null(fit) || isTRUE(fit$nevent == 0)) -Inf else {
                           target_coef <- paste0("cut_design", numcut + 1)
                           if (!(target_coef %in% names(fit$coefficients)) || is.na(fit$coefficients[target_coef])) -Inf else exp(fit$coefficients[target_coef])
                         }
                       },
                       "p_value" = {
                         fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
                         if (is.null(fit) || is.null(fit$loglik) || isTRUE(fit$nevent == 0)) -Inf else {
                           if (!is.na(loglik0)) {
                             pval <- stats::pchisq(-2 * (loglik0 - fit$loglik[2]), numcut, lower.tail = FALSE)
                             1 - pval
                           } else {
                             sfit <- tryCatch(summary(fit), error = function(e) NULL)
                             if (is.null(sfit) || is.null(sfit$logtest)) -Inf else 1 - sfit$logtest["pvalue"]
                           }
                         }
                       },
                       "loglik" = {
                         beta_init <- params[(numcut + 1):length(params)]

                         # DEFENSIVE ENGINE BACKSTOP: If running an unadjusted model,
                         # bypass manual init arrays to let coxph resolve loglik natively
                         if (length(beta_init) == 0 || is.null(confound) || ncol(confound) == 0) {
                           fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
                           if (is.null(fit) || is.null(fit$loglik) || isTRUE(fit$nevent == 0)) -Inf else fit$loglik[2]
                         } else {
                           fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit, init = beta_init, iter.max = 0), error = function(e) NULL)
                           if (is.null(fit) || is.null(fit$loglik) || isTRUE(fit$nevent == 0)) -Inf else fit$loglik[2]
                         }
                       },
                       -Inf
  )

  if (!is.null(cache)) cache[[key]] <- stat_value
  return(stat_value)
}

#' Internal helper: Wrapper for the rgenoud Genetic Algorithm
#'
#' @description
#' Sets up and runs `rgenoud::genoud` to find optimal cut-points, defining
#' search domains, caching environments, and initial starting values.
#'
#' @param target Continuous predictor vector.
#' @param numcut Number of cut-points.
#' @param time Survival time vector.
#' @param censor Survival event vector.
#' @param confound Optional covariate data frame.
#' @param nmin Minimum observations per group.
#' @param criterion Optimisation criterion.
#' @param max.generations Max generations (native rgenoud arg).
#' @param pop.size Population size (native rgenoud arg).
#' @param gap Minimum distance between cut-points.
#' @param print.level Console output level.
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
#' @srrstats {G2.0} `requireNamespace("rgenoud")` check.
#' @srrstats {G5.8} Edge cases (zero gap, infinite range) return `NULL`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom stats quantile as.formula na.omit
#' @importFrom survival coxph Surv
#' @noRd
.run_genetic_search <- function(target, numcut, time, censor, confound, nmin,
                                criterion, max.generations = 100, pop.size = 100,
                                gap = NULL, print.level = 0, ...) {
  if (!requireNamespace("rgenoud", quietly = TRUE)) {
    stop("The 'rgenoud' package is required. Please install it.", call. = FALSE)
  }

  num_confound_vars <- if (is.null(confound)) 0 else ncol(confound)
  optimising_betas <- (criterion == "loglik")
  nvars <- if (optimising_betas) numcut * 2 + num_confound_vars else numcut

  domain_cuts <- matrix(rep(range(target, na.rm = TRUE), numcut), ncol = 2, byrow = TRUE)
  if (any(is.infinite(domain_cuts))) return(NULL)
  domain <- if (optimising_betas) rbind(domain_cuts, matrix(rep(c(-5, 5), nvars - numcut), ncol = 2, byrow = TRUE)) else domain_cuts

  initial_cuts <- stats::quantile(target, probs = seq(0, 1, length.out = numcut + 2), na.rm = TRUE)[2:(numcut + 1)]
  initial_values <- if (optimising_betas) c(initial_cuts, rep(0, nvars - numcut)) else initial_cuts

  if (is.null(gap)) {
    gap <- stats::quantile(sort(diff(sort(unique(stats::na.omit(target))))), probs = 0.5, na.rm = TRUE)
    if (is.na(gap) || gap == 0) gap <- 1e-4
  }

  loglik0 <- NA_real_
  if (criterion == "p_value") {
    null_formula_str <- "survival::Surv(time, censor) ~ 1"
    data_for_null_fit <- data.frame(time = time, censor = censor)
    if (num_confound_vars > 0) {
      colnames(confound) <- make.names(colnames(confound))
      data_for_null_fit <- cbind(data_for_null_fit, confound)
      null_formula_str <- paste(null_formula_str, "+", paste(colnames(confound), collapse = " + "))
    }
    null_fit <- tryCatch(survival::coxph(stats::as.formula(null_formula_str), data = data_for_null_fit), error = function(e) NULL)
    if (!is.null(null_fit)) loglik0 <- null_fit$loglik[2]
  }

  eval_cache <- new.env(hash = TRUE, parent = emptyenv())
  base_df <- data.frame(time = time, censor = censor)
  formula_str <- "survival::Surv(time, censor) ~ cut_design"
  if (num_confound_vars > 0) {
    colnames(confound) <- make.names(colnames(confound))
    base_df <- cbind(base_df, confound)
    formula_str <- paste(formula_str, "+", paste(colnames(confound), collapse = " + "))
  }

  optim_result <- tryCatch(rgenoud::genoud(
    fn = .obj, nvars = nvars, max = TRUE, pop.size = pop.size,
    max.generations = max.generations, wait.generations = 10, hard.generation.limit = TRUE,
    starting.values = initial_values, Domains = domain, print.level = print.level,
    time = time, censor = censor, target = target, confound = confound,
    numcut = numcut, gap = gap, nmin = nmin, criterion = criterion,
    loglik0 = loglik0, cache = eval_cache, base_df = base_df,
    precompiled_formula = stats::as.formula(formula_str)
  ), error = function(e) NULL)

  return(optim_result)
}

#' Internal helper: Genetic Model Selection for finding `max_cuts`
#'
#' @description
#' Evaluates models from 1 to `max_cuts` using the `rgenoud` evolutionary engine,
#' computing Information Criteria (AIC/BIC) to penalize complexity.
#'
#' @param userdata Cleaned survival data frame.
#' @param max_cuts Maximum number of cut-points to test.
#' @param nmin Minimum observations per group.
#' @param criterion IC to calculate (AIC, AICc, BIC).
#' @param covariates Optional covariates.
#' @param max.generations Iterations for rgenoud.
#' @param pop.size Population size for rgenoud.
#' @param ... Unused arguments safely passed down.
#'
#' @return A data.frame of model selection results.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Helper for genetic-algorithm-based model selection.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @importFrom cli cli_inform cli_alert_info
#' @noRd
.genetic_search_num <- function(userdata, max_cuts, nmin, criterion, covariates, max.generations, pop.size, ...) {
  n <- nrow(userdata)
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""
  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)

  ic0 <- tryCatch({
    fit0 <- survival::coxph(stats::as.formula(base_formula_str), data = userdata)
    .calc_ic(fit0, k = 1 + length(covariates), n = n, criterion = criterion)
  }, error = function(e) {
    if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("Could not calculate IC for base model (0 cuts):", e$message))
    NA_real_
  })

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    if (requireNamespace("cli", quietly = TRUE)) cli::cli_alert_info(paste("Running genetic algorithm for", k_cuts, "cut-point(s)..."))

    ga_result <- tryCatch({
      .run_genetic_search(
        target = userdata$factor, numcut = k_cuts, time = userdata$time, censor = userdata$event,
        confound = if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL,
        nmin = nmin, criterion = "loglik", max.generations = max.generations, pop.size = pop.size, ...
      )
    }, error = function(e) {
      if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("Genetic algorithm failed for", k_cuts, "cut(s):", e$message))
      NULL
    })

    if (!is.null(ga_result) && is.finite(ga_result$value) && ga_result$value > -.Machine$double.xmax) {
      ic_val <- .calc_ic(list(loglik = c(NA, ga_result$value)), k = k_cuts + length(covariates), n = n, criterion = criterion)
      results <- rbind(results, data.frame(num_cuts = k_cuts, IC = ic_val, cuts = I(list(sort(ga_result$par[1:k_cuts])))))
    } else {
      if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("No valid cut-points found for", k_cuts, "cut(s) due to genetic algorithm failure or constraints."))
      results <- rbind(results, data.frame(num_cuts = k_cuts, IC = NA_real_, cuts = I(list(NULL))))
    }
  }
  names(results)[2] <- criterion
  return(results)
}

#' Internal helper: Calculate Information Criterion (AIC, AICc, or BIC)
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
#' @srrstats {RE4.11} Implements AIC, AICc, BIC.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.calc_ic <- function(model, k, n, criterion) {
  if (is.null(model) || !is.list(model) || is.null(model$loglik)) {
    return(NA_real_)
  }
  logL <- model$loglik[2]
  if (is.na(logL)) {
    return(NA_real_)
  }
  if (criterion == "BIC" && (is.na(n) || n <= 0)) {
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
  } else {
    return(-2 * logL + 2 * k)
  }
}
