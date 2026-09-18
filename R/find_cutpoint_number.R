# ===================================================================
#' Find Optimal Number of Cut-points for Survival Data
# ===================================================================
#' @description
#' Finds optimal cut-point number (0 to `max_cuts`) for a Cox model
#' by comparing AIC, AICc, or BIC. Features hardware-accelerated grouping
#' iterations via Rcpp compilation hooks and robust UX constraint warnings.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Computes Goodness-of-fit (AIC, AICc, BIC) for selection.
#' @srrstats {G2.13} `cli_abort()` checks inputs pre-processing.
#' @srrstats {G2.14} `NA` results handled via `na_result` helper.
#' @srrstats {RE6.0} `plot()` method provided.
#' @srrstats {RE6.2} `plot()` visualises model selection metric vs cuts.
#' @srrstats {G1.0} References provided for AIC, AICc, BIC.
#' @srrstats {G1.3} Systematic grid search (max_cuts <= 2) and
#' `rgenoud` global optimisation documented.
#' @srrstats {G1.5} Compared with `cutpointr`/`survminer` in vignette.
#' @srrstats {G1.6} Numerical stability via `survival::coxph`
#' and `rgenoud`; edge cases return `NA`.
#' @srrstats {G2.3b} Uses `as.formula()`/`data[...]`; no unsafe eval.
#' @srrstats {G2.4} `NA` removed via `stats::na.omit()`.
#' @srrstats {G2.4b} Explicit conversion checks prevent factors or character strings from causing silent optimisation failures.
#' @srrstats {G5.2} `optimal_num_cuts` is `NA` when no valid solution found.
#' @srrstats {G2.6} Input validation via direct checks.
#' @srrstats {G2.8} Informative errors via `cli::cli_abort()`.
#' @srrstats {G5.2} Warnings via `cli::cli_alert_warning()`.
#' @srrstats {G5.2} Graceful degradation via `na_result()` for
#' empty data or model failures.
#' @srrstats {G2.14c} `NA` propagation controlled.
#' @srrstats {G1.4} All parameters/return values documented.
#' @srrstats {RE4.17} `print()` method provided.
#' @srrstats {RE4.18} `summary()` method provided.
#'
#' @srrstats {RE1.1} Assumes PH; check `summary()` for `cox.zph`.
#' @srrstats {RE1.3} PH diagnostics via `cox.zph` in `summary()`.
#' @srrstats {RE1.3a} `summary()` includes PH test results.
#' @srrstats {RE2.0} Estimates/SEs from `coxph` in `summary()`.
#' @srrstats {RE2.1} Missing values handled via explicit `na.omit`.
#' @srrstats {RE3.0} `tryCatch` checks model convergence; failures return `NA`.
#' @srrstats {RE2.4a} Checks for collinearity via model fitting constraints.
#' @srrstats {RE2.4b} Checks for insufficient data/constant predictor.
#' @srrstats {RE4.2} Model selection via AIC, AICc, or BIC.
#' @srrstats {RE5.0} Model averaging not implemented.
#' @srrstats {RE6.3} No diagnostic plots; use `cox.zph`.
#' @srrstats {RE7.0a} `na.omit()` removes missing data.
#'
#' @details
#' `method = "systematic"`: grid search respecting `nmin`.
#' `method = "genetic"`: `rgenoud` global optimisation.
#' Systematic search is slow for `max_cuts > 2`; use `genetic`.
#' Core vector partitions are calculated in compiled C++ via `Rcpp` for optimal performance.
#'
#' @param data Input data frame.
#' @param predictor Continuous predictor variable name (character).
#' @param outcome_time Time-to-event variable name (character).
#' @param outcome_event Event indicator name (0/1) (character).
#' @param method `"systematic"` (max_cuts <= 2) or `"genetic"`.
#' @param criterion `"AIC"`, `"AICc"` or `"BIC"`.
#' @param covariates Character vector of covariate names (optional).
#' @param max_cuts Max number of cut-points to test (non-negative int).
#' @param nmin Min. group size (count or proportion).
#' @param seed Integer or `NULL`; random seed for `rgenoud`.
#' @param max.generations Integer; generations for `rgenoud`. If `NULL`, dynamically scales.
#' @param pop.size Integer; population size for `rgenoud`. If `NULL`, dynamically scales.
#' @param use_cpp Logical. Automatically checks and calls compiled C++ routines via `Rcpp`. Default is `TRUE`.
#' @param x An object from [find_cutpoint_number()].
#' @param object An object from [find_cutpoint_number()].
#' @param y Unused.
#' @param ... Additional arguments passed to `rgenoud`.
#'
#' @return An S3 object (`find_cutpoint_number_result`) with
#' `results`, `parameters`, `userdata`, `optimal_num_cuts`,
#' `optimal_cuts`, and `candidate_cuts`.
#'
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom stats na.omit as.formula aggregate quantile symnum
#' @importFrom survival coxph Surv survfit
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success
#' @importFrom cli cli_inform cli_abort cli_bullets cli_h2 cli_alert_warning
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#' @importFrom ggplot2 theme_minimal scale_x_continuous element_text theme
#' @importFrom tools toTitleCase
#' @export
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'
#'   set.seed(42)
#'   n <- 100
#'   biomarker <- rnorm(n, mean = 6, sd = 1.2)
#'   sim_data <- data.frame(
#'     time      = rexp(n, rate = 0.05 + 0.03 * (biomarker > 6)),
#'     event     = rbinom(n, 1, 0.8),
#'     biomarker = biomarker
#'   )
#'
#'   num_fit <- find_cutpoint_number(
#'     data          = sim_data,
#'     predictor     = "biomarker",
#'     outcome_time  = "time",
#'     outcome_event = "event",
#'     max_cuts      = 1,
#'     method        = "systematic",
#'     criterion     = "BIC",
#'     nmin          = 0.2,
#'     quiet         = TRUE
#'   )
#'   summary(num_fit)
#' }
find_cutpoint_number <- function(data, predictor,
                                 outcome_time, outcome_event,
                                 method = "systematic", criterion = "BIC",
                                 covariates = NULL, max_cuts = 2,
                                 nmin = 0.1, seed = NULL,
                                 max.generations = NULL, pop.size = NULL,
                                 use_cpp = TRUE, ...) {
  if (!is.numeric(max_cuts) || max_cuts < 0 || max_cuts != round(max_cuts)) {
    cli::cli_abort("max_cuts must be a non-negative integer.")
  }
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("BIC", "AIC", "AICc"))
  if (is.null(predictor)) {
    cli::cli_abort("A 'predictor' variable must be specified.")
  }
  if (is.null(outcome_time) || is.null(outcome_event)) {
    cli::cli_abort("'outcome_time' and 'outcome_event' are required.")
  }

  if (method == "genetic" && !requireNamespace("rgenoud", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "'genetic' method requires the 'rgenoud' package.",
        "i1" = "Install with: install.packages(\"rgenoud\")",
        "i2" = "Or, use `method = \"systematic\"` (for max_cuts <= 2)."
      )
    )
  }

  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  if (!all(required_vars %in% names(data))) {
    missing_cols <- required_vars[!required_vars %in% names(data)]
    cli::cli_abort(
      "Missing columns: {paste(missing_cols, collapse = ', ')}"
    )
  }

  original_predictor_name <- predictor

  userdata <- data[, required_vars, drop = FALSE]
  userdata <- stats::na.omit(userdata)

  if (is.factor(userdata[[predictor]]) || is.character(userdata[[predictor]])) {
    cli::cli_abort(c(
      "x" = "Unordered factors or characters are not allowed as continuous optimisation inputs.",
      "i" = "Variable '{predictor}' must be passed as a continuous numeric vector."
    ))
  }

  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"

  n <- nrow(userdata)

  .validate_event_column(userdata$event, outcome_event)

  if (use_cpp && !exists("cpp_get_group_assignments", mode = "function")) {
    cli::cli_alert_warning("Compiled C++ binary not loaded. Falling back gracefully to native R processing.")
    use_cpp <- FALSE
  }

  if (nmin < 1 && nmin > 0) {
    nmin_abs <- floor(nmin * n)
    cli::cli_alert_info("nmin {nmin} is a proportion. Min. group size set to {nmin_abs}.")
  } else if (nmin >= 1) {
    nmin_abs <- as.integer(nmin)
  } else {
    cli::cli_abort("'nmin' must be a positive number.")
  }

  # ✅ FIXED: Keep ... in the function signature if needed, but drop it from the list constructor below
  na_result <- function(userdata, ...) {
    output <- list(
      results = data.frame(),
      parameters = list(
        method = method,
        criterion = criterion,
        analysis_type = "survival",
        predictor = original_predictor_name,
        outcome_time = outcome_time,
        outcome_event = outcome_event,
        covariates = covariates,
        max_cuts = max_cuts,
        nmin = nmin_abs,
        max.generations = max.generations,
        pop.size = pop.size,
        use_cpp = use_cpp
      ),
      userdata = userdata,
      optimal_num_cuts = NA,
      optimal_cuts = NA,
      candidate_cuts = NULL
    )
    class(output) <- "find_cutpoint_number_result"
    cli::cli_inform("No valid results found with given parameters.")
    return(output)
  }

  if (n == 0) {
    cli::cli_inform("No complete cases found after removing NAs.")
    return(na_result(userdata))
  }

  if (n < nmin_abs * (max_cuts + 1)) {
    cli::cli_alert_warning(
      "Sample size ({n}) is mathematically insufficient to support {max_cuts} cuts with a minimum group size constraint of {nmin_abs} patients ({max_cuts + 1} groups required)."
    )
    return(na_result(userdata))
  }

  # Build the base Regularised Grid Space lattice to gauge unique data supply
  grid_probs_vec <- seq(0.01, 0.99, by = 0.01)
  permissible_cuts <- sort(unique(stats::quantile(userdata$factor, probs = grid_probs_vec, na.rm = TRUE)))

  # ===================================================================
  # UX COHORT HEADROOM CHECK (Dynamic Lot Defense)
  # ===================================================================
  unique_values_needed <- (max_cuts + 1) * nmin_abs
  actual_available_nodes <- length(permissible_cuts)

  if (n < unique_values_needed || actual_available_nodes < max_cuts) {
    cli::cli_abort(c(
      "x" = "Insufficient data density for {max_cuts} cut-point(s) at nmin = {nmin_abs}.",
      "i" = "Model requires at least {max_cuts} unique partition boundaries, but only {actual_available_nodes} survived duplicate filtering.",
      "*" = "Action: Reduce 'max_cuts' or lower your 'nmin' threshold."
    ))
  }

  if (length(unique(userdata$factor)) <= max_cuts) {
    cli::cli_inform(paste(
      "Predictor has too few unique values",
      "({length(unique(userdata$factor))}) for max_cuts ({max_cuts})."
    ))
    return(na_result(userdata))
  }

  cli::cli_alert_info("Finding optimal cut number: method = {.strong {method}}")

  if (method == "systematic") {
    params <- list(
      userdata = userdata, max_cuts = max_cuts, nmin = nmin_abs,
      criterion = criterion, covariates = covariates,
      max.generations = max.generations, pop.size = pop.size,
      use_cpp = use_cpp
    )
    extra_args <- list(...)
    params <- c(params, extra_args)

    search_payload <- do.call(.systematic_search_num, params)
    results <- search_payload$selection_table
    permissible_cuts <- search_payload$candidate_cuts
  } else { # genetic complexity loop sweep
    results_list <- list()

    results_list[[1]] <- data.frame(
      num_cuts = 0,
      BIC = .get_model_ic_num(userdata, rep(1, n), 0, n, criterion, covariates),
      cuts = I(list(NA))
    )
    names(results_list[[1]])[2] <- criterion

    confound_df <- if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL

    for (k in 1:max_cuts) {
      cli::cli_inform("Running discrete genetic algorithm for {k} cut-point(s)...")

      # --- SYSTEM SYNCHRONISED MEMORY-SAFE TUNING ---
      current_pop <- pop.size
      current_gen <- max.generations
      if (is.null(pop.size) || pop.size == 100) {
        current_pop <- switch(as.character(k),
          "1" = 30,
          "2" = 60,
          "3" = 100,
          "4" = 120,
          "5" = 150,
          "6" = 180,
          "7" = 200,
          250
        )
      }
      if (is.null(max.generations) || max.generations == 100) {
        current_gen <- switch(as.character(k),
          "1" = 30,
          "2" = 40,
          "3" = 50,
          "4" = 55,
          "5" = 60,
          "6" = 65,
          "7" = 70,
          80
        )
      }

      params_k <- list(
        target = userdata$factor,
        numcut = k,
        time = userdata$time,
        censor = userdata$event,
        confound = confound_df,
        nmin = nmin_abs,
        criterion = if (criterion %in% c("BIC", "AIC", "AICc")) "p_value" else criterion,
        max.generations = current_gen,
        pop.size = current_pop,
        use_cpp = use_cpp,
        print.level = 0,
        grid_probs = grid_probs_vec
      )
      extra_args <- list(...)

      if (!("boundary.enforcement" %in% names(extra_args))) {
        extra_args$boundary.enforcement <- 2
      }
      params_k <- c(params_k, extra_args)

      raw_fit <- do.call(".run_genetic_search", params_k)

      if (is.null(raw_fit) || !is.finite(raw_fit$value) || raw_fit$value <= -.Machine$double.xmax) {
        fit_k <- NULL
        # DYNAMIC RECOGNITION ARREST: Throw warning for dropped candidate configurations
        cli::cli_alert_warning("Model with {k} cut-point(s) collapsed: Subgroups violated the minimum 'nmin' constraint of {nmin_abs} subjects.")
      } else {
        fit_k <- data.frame(
          num_cuts = k,
          placeholder = .get_model_ic_num(userdata, cut(userdata$factor, breaks = c(-Inf, sort(raw_fit$par[1:k]), Inf), labels = FALSE), k, n, criterion, covariates),
          cuts = I(list(sort(raw_fit$par[1:k])))
        )
        names(fit_k)[2] <- criterion
      }

      if (!is.null(fit_k)) {
        results_list[[k + 1]] <- fit_k
      }
    }
    results <- do.call(rbind, results_list)
  }

  if (is.null(results) || !is.data.frame(results) || nrow(results) == 0) {
    cli::cli_inform("Search algorithm failed to produce results.")
    return(na_result(userdata))
  }

  if (anyNA(results[[criterion]])) {
    cli::cli_alert_warning("All tested model cut-points violated localised subgroup size constraints during runtime search iterations.")
  }

  min_ic <- min(results[[criterion]], na.rm = TRUE)
  delta_col_name <- paste0("Delta_", criterion)
  weight_col_name <- paste0(criterion, "_Weight")

  results[[delta_col_name]] <- results[[criterion]] - min_ic
  exp_delta <- exp(-0.5 * results[[delta_col_name]])
  results[[weight_col_name]] <- exp_delta / sum(exp_delta, na.rm = TRUE)

  results$Evidence <- vapply(results[[delta_col_name]], function(d) {
    if (is.na(d)) {
      NA_character_
    } else if (d <= 2) {
      "Substantial"
    } else if (d <= 7) {
      "Moderate"
    } else {
      "Minimal"
    }
  }, FUN.VALUE = character(1))

  output <- list(
    results = results,
    parameters = list(
      method = method,
      criterion = criterion,
      analysis_type = "survival",
      predictor = original_predictor_name,
      outcome_time = outcome_time,
      outcome_event = outcome_event,
      covariates = covariates,
      max_cuts = max_cuts,
      nmin = nmin_abs,
      max.generations = max.generations,
      pop.size = pop.size,
      use_cpp = use_cpp
    ),
    userdata = userdata,
    candidate_cuts = permissible_cuts
  )

  finite_ic <- results[[criterion]][is.finite(results[[criterion]])]
  if (length(finite_ic) > 0) {
    min_ic_idx <- which.min(results[[criterion]])
    output$optimal_num_cuts <- results$num_cuts[min_ic_idx]
    output$optimal_cuts <- results$cuts[[min_ic_idx]]
  } else {
    output$optimal_num_cuts <- NA
    output$optimal_cuts <- NA
  }

  class(output) <- "find_cutpoint_number_result"
  return(output)
}

#' Internal helper: Compute IC from standard factors for number selection
#'
#' @description
#' Computes the Information Criterion dynamically for the selection engine
#' based on the evaluated grouping factors, defensively capturing model singularities.
#'
#' @inheritParams find_cutpoint_number
#' @param factor_status The cut survival factor.
#' @param k_cuts Number of cuts evaluated.
#' @param n Sample size.
#' @param cov_part Character vector of covariate column names.
#'
#' @return Single numeric IC value (or `Inf` on failure).
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE3.1} Singularity and non-convergence exceptions in the variance-covariance matrix are gracefully intercepted via tryCatch blocks.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @noRd
.get_model_ic_num <- function(userdata, factor_status, k_cuts, n,
                              criterion, cov_part) {
  num_cov <- length(cov_part[cov_part != ""])
  cov_string <- if (num_cov > 0) paste(" +", paste(cov_part, collapse = " + ")) else ""
  formula_str <- paste("survival::Surv(time, event) ~ factor_status", cov_string)

  fit <- tryCatch(
    survival::coxph(as.formula(formula_str), data = userdata),
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(fit)) {
    return(Inf)
  }

  k_params <- k_cuts + num_cov
  .calc_ic(fit, k_params, n, criterion)
}

# print.find_cutpoint_number_result(), summary.find_cutpoint_number_result(),
# and plot.find_cutpoint_number_result() were previously duplicated here with
# a different (older, less complete) implementation than the versions in
# find_cutpoint_number_methods.R - missing the Schoenfeld/proportional-hazards
# section, the zero-event safeguard, the all.x = TRUE merge fix, and using
# ggplot2::theme_minimal() instead of the package's theme_optsurv(). The
# copies in find_cutpoint_number_methods.R are confirmed to be the versions
# actually in effect (consistent with observed console/plot output across
# this package's case studies and vignettes) and are now the single source.
