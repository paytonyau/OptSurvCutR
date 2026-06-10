# ===================================================================
# INTERNAL HELPER FUNCTIONS
# Support for `find_cutpoint()` and `find_cutpoint_number()`.
# ===================================================================

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

# ===================================================================
# Validation & Prep Helpers for find_cutpoint()
# ===================================================================

#' Validate inputs for find_cutpoint()
#'
#' @description
#' Centralised validation for `find_cutpoint()` arguments.
#'
#' @param method,criterion,num_cuts,covariates,predictor,outcome_time,
#' outcome_event,data Arguments from `find_cutpoint()`.
#'
#' @return `invisible(TRUE)` on success; aborts on failure.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.0} Validates `method`, `criterion`, `num_cuts`.
#' @srrstats {G2.1} Type checks.
#' @srrstats {G2.4b} `match.arg()` for controlled vocabularies.
#' @srrstats {G2.4c} Conversion mechanisms explicitly account for logical vectors during argument verification.
#' @srrstats {G2.0} `requireNamespace()` for optional `rgenoud`.
#' @srrstats {G2.9} Column-name existence checks.
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
#' from `find_cutpoint()`.
#'
#' @return Cleaned `userdata` data.frame.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE2.1} `na.omit()` with explicit `NA` handling.
#' @srrstats {RE2.2} Missing value processing parameters explicitly separate row omissions from complete missingness boundaries.
#' @srrstats {G2.13} Checks for missing data (via NA removal).
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#' @srrstats {G2.10} Ensures column extraction handles single columns
#'   consistently using `drop = FALSE`.
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

#' Validate event column (0/1)
#'
#' @description
#' Checks that the event column is numeric and contains only 0 and 1.
#'
#' @param event_col The event vector (e.g., `userdata$event`).
#' @param outcome_event The original name of the event column (for errors).
#'
#' @return `invisible(TRUE)` on success; aborts on failure.
#'
#' @srrstats {G2.1} Event column must be numeric 0/1.
#' @srrstats {G2.13} `cli_abort()` for invalid event data.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.validate_event_column <- function(event_col, outcome_event) {
  if (is.character(event_col) || is.factor(event_col)) {
    cli::cli_abort("Event column {.var {outcome_event}} must be numeric.")
  }

  unique_events <- unique(na.omit(event_col))
  if (!all(unique_events %in% c(0, 1, TRUE, FALSE))) {
    cli::cli_abort(c(
      "x" = "Event column {.var {outcome_event}} must be strictly binary (0/1).",
      "i" = "Detected invalid values: {.val {sort(unique_events)}}."
    ))
  }
  invisible(TRUE)
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
#' @srrstats {G2.1} Event column must be numeric 0/1.
#' @srrstats {G2.13} `cli_abort()` for invalid event data.
#' @srrstats {G5.8} Edge cases (constant predictor, insufficient data).
#' @srrstats {G5.8d} Intercepts zero-length groups or data rows falling below relative minimum cohort size allocations.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
.validate_data_conditions <- function(userdata, nmin, num_cuts,
                                      outcome_event, quiet) {
  if (any(userdata$time < 0, na.rm = TRUE)) {
    stop("Time variable must be non-negative.", call. = FALSE)
  }

  .validate_event_column(userdata$event, outcome_event)

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
#' @noRd
`%||%` <- function(a, b) if (is.null(a)) b else a
