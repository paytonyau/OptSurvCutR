# ===================================================================
# INTERNAL ENGINE: SYSTEMATIC SEARCH
# Handles exhaustive grid searches for find_cutpoint and find_cutpoint_number
# ===================================================================

#' Internal helper: Systematic Grid Search
#'
#' @description
#' Implements an exhaustive grid search to evaluate all possible thresholds
#' for 1 or 2 cut-points, respecting the minimum group size constraints.
#'
#' @param userdata Cleaned survival data frame.
#' @param num_cuts Number of cut-points to evaluate (1 or 2).
#' @param criterion Statistic to optimise: "logrank", "hazard_ratio", or "p_value".
#' @param covariates Optional vector of covariate names.
#' @param nmin Absolute minimum number of observations per group.
#' @param predictor_name Original name of the predictor (for messaging).
#' @param quiet Logical to suppress console output.
#' @param ... Additional unused arguments (absorbed safely).
#'
#' @return A list containing `optimal_cuts`, `optimal_stat`, and `all_stats`.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE1.0} Implements systematic grid search for 1–2 cut-points.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom cli cli_alert_info cli_alert_warning cli_inform cli_alert_success
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @noRd
.systematic_search <- function(userdata, num_cuts, criterion,
                               covariates, nmin, predictor_name,
                               quiet, ...) {
  if (!quiet) cli::cli_alert_info("Running systematic search for {num_cuts} cut-point(s)...")
  userdata <- userdata[order(userdata$factor), ]

  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

  fit_null_model <- NULL
  if (criterion == "p_value") {
    null_formula_str <- if (!is.null(covariates)) {
      paste("Surv(time, event) ~", paste(covariates, collapse = " + "))
    } else {
      "Surv(time, event) ~ 1"
    }

    if (length(unique(userdata$time)) <= 1) {
      fit_null_model <- NULL
    } else {
      fit_null_model <- tryCatch(
        survival::coxph(stats::as.formula(null_formula_str), data = userdata),
        error = function(e) NULL
      )
    }

    if (is.null(fit_null_model)) {
      if (!quiet) cli::cli_alert_warning("Could not fit null model for p-value. Aborting.")
      return(list(
        optimal_cuts = rep(NA_real_, num_cuts), optimal_stat = NA_real_,
        all_stats = NULL, parameters = list(method = "systematic")
      ))
    }
  }

  direction <- if (criterion == "p_value") "min" else "max"
  best_stat <- if (direction == "min") Inf else -Inf
  best_cut_val <- rep(NA_real_, num_cuts)
  all_stats_df <- NULL

  if (num_cuts == 1) {
    search_grid <- unique(userdata$factor[nmin:(nrow(userdata) - nmin)])
    if (length(search_grid) == 0) {
      return(list(optimal_cuts = NA_real_, optimal_stat = NA_real_, all_stats = NULL, parameters = list(method = "systematic")))
    }

    stats_per_cut <- vapply(search_grid, .get_stat,
      num_cuts = 1, data_in = userdata,
      criterion = criterion, cov_formula = cov_part,
      nmin = nmin, fit_null = fit_null_model,
      FUN.VALUE = numeric(1)
    )

    if (all(is.na(stats_per_cut))) {
      return(list(optimal_cuts = NA_real_, optimal_stat = NA_real_, all_stats = NULL, parameters = list(method = "systematic")))
    }

    best_idx <- if (direction == "min") which.min(stats_per_cut) else which.max(stats_per_cut)
    best_cut_val <- search_grid[best_idx]
    best_stat <- stats_per_cut[best_idx]
    all_stats_df <- data.frame(cut1 = search_grid, stat = stats_per_cut)
  } else { # num_cuts == 2
    if (!quiet) cli::cli_alert_info("Searching for 2 cuts is slow...")
    foreach::registerDoSEQ()

    possible_c1_indices <- nmin:(nrow(userdata) - (2 * nmin))
    grid1_values <- unique(userdata$factor[possible_c1_indices])

    if (length(grid1_values) == 0) {
      return(list(optimal_cuts = c(NA_real_, NA_real_), optimal_stat = NA_real_, all_stats = NULL, parameters = list(method = "systematic")))
    }

    results_list <- foreach::foreach(c1 = grid1_values, .combine = "rbind", .export = c(".get_stat")) %do% {
      best_local_stat <- if (direction == "min") Inf else -Inf
      best_local_c2 <- NA_real_
      c1_max_index <- max(which(userdata$factor == c1))

      start_index_c2 <- c1_max_index + nmin
      end_index_c2 <- nrow(userdata) - nmin
      if (start_index_c2 > end_index_c2) {
        return(NULL)
      }

      grid2_values <- unique(userdata$factor[start_index_c2:end_index_c2])
      if (length(grid2_values) == 0) {
        return(NULL)
      }

      for (c2 in grid2_values) {
        stat <- .get_stat(c(c1, c2), 2, userdata, criterion, cov_part, nmin, fit_null = fit_null_model)
        if (is.na(stat)) next
        is_better <- if (direction == "min") (stat < best_local_stat) else (stat > best_local_stat)
        if (is_better && !is.infinite(stat)) {
          best_local_stat <- stat
          best_local_c2 <- c2
        }
      }
      if (is.na(best_local_c2)) {
        return(NULL)
      }
      data.frame(stat = best_local_stat, c1 = c1, c2 = c2)
    }

    if (is.null(results_list) || nrow(results_list) == 0) {
      return(list(optimal_cuts = c(NA_real_, NA_real_), optimal_stat = NA_real_, all_stats = NULL, parameters = list(method = "systematic")))
    }

    best_idx <- if (direction == "min") which.min(results_list$stat) else which.max(results_list$stat)
    best_stat <- results_list$stat[best_idx]
    best_cut_val <- c(results_list$c1[best_idx], results_list$c2[best_idx])
  }

  if (!quiet) cli::cli_alert_success("Systematic search complete.")
  return(list(
    optimal_cuts = best_cut_val,
    optimal_stat = best_stat,
    all_stats = all_stats_df,
    parameters = list(method = "systematic") # Satisfies internal tests
  ))
}

#' Internal helper: Calculate Survival Statistic
#'
#' @description
#' Computes the target statistic (Log-Rank, Hazard Ratio, or P-value) for a
#' specific set of proposed cut-points.
#'
#' @param cuts Vector of proposed cut-points.
#' @param num_cuts Number of cut-points.
#' @param data_in Sub-setting of the survival data.
#' @param criterion Statistic to compute.
#' @param cov_formula Covariate string for the formula.
#' @param nmin Minimum observations per group.
#' @param fit_null Pre-computed null model for LRT p-value calculation.
#'
#' @return Single numeric statistic value or NA.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE1.0} Uses `survdiff`/`coxph` for log-rank, HR, p-value.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival survdiff coxph Surv
#' @importFrom stats as.formula pchisq
#' @noRd
.get_stat <- function(cuts, num_cuts, data_in, criterion, cov_formula, nmin, fit_null = NULL) {
  if (length(unique(data_in$time)) <= 1) {
    return(NA)
  }
  sorted_cuts <- sort(cuts)
  data_in$group <- as.factor(findInterval(data_in$factor, sorted_cuts, left.open = TRUE) + 1L)

  if (any(table(data_in$group) < nmin) || nlevels(data_in$group) != (num_cuts + 1)) {
    return(NA)
  }

  formula_str <- paste("Surv(time, event) ~ group", cov_formula)

  if (criterion == "logrank") {
    if (cov_formula == "") {
      fit <- tryCatch(survival::survdiff(stats::as.formula(formula_str), data = data_in), error = function(e) NULL)
      if (is.null(fit)) {
        return(NA)
      }
      return(fit$chisq)
    } else {
      fit <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = data_in), error = function(e) NULL)
      if (is.null(fit) || is.null(fit$score)) {
        return(NA)
      }
      return(fit$score)
    }
  } else {
    fit <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = data_in), error = function(e) NULL)
    if (is.null(fit)) {
      return(NA)
    }

    if (criterion == "hazard_ratio") {
      if (is.null(fit$coefficients)) {
        return(NA)
      }
      coef_name <- paste0("group", num_cuts + 1)
      if (!(coef_name %in% names(fit$coefficients)) || is.na(fit$coefficients[coef_name])) {
        return(-Inf)
      }
      return(exp(fit$coefficients[coef_name]))
    } else if (criterion == "p_value") {
      if (is.null(fit_null) || is.null(fit$loglik)) {
        return(NA)
      }
      lrt_stat <- 2 * (fit$loglik[2] - fit_null$loglik[2])
      return(stats::pchisq(lrt_stat, df = num_cuts, lower.tail = FALSE))
    }
  }
}

#' Internal helper: Systematic Model Selection for finding `max_cuts`
#'
#' @description
#' Evaluates models from 1 to `max_cuts` using the systematic exhaustive
#' engine, computing Information Criteria (AIC/BIC) for model selection.
#'
#' @param userdata Cleaned survival data frame.
#' @param max_cuts Maximum number of cut-points to test (<= 2).
#' @param nmin Minimum observations per group.
#' @param criterion IC to calculate (AIC, AICc, BIC).
#' @param covariates Optional covariates.
#' @param ... Unused arguments passed down safely.
#'
#' @return A data.frame of model selection results.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Helper for systematic model selection (grid search).
#' @srrstats {RE1.0} Fits Cox models via `survival::coxph`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @importFrom cli cli_inform cli_alert_info
#' @noRd
.systematic_search_num <- function(userdata, max_cuts, nmin, criterion, covariates, ...) {
  if (max_cuts > 2) stop("systematic method is only implemented for max_cuts <= 2.")
  n <- nrow(userdata)
  userdata <- userdata[order(userdata$factor), ]
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)
  ic0 <- tryCatch(
    {
      fit0 <- survival::coxph(stats::as.formula(base_formula_str), data = userdata)
      .calc_ic(fit0, k = 1 + length(covariates), n = n, criterion = criterion)
    },
    error = function(e) {
      if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("Could not calculate IC for base model (0 cuts):", e$message))
      NA_real_
    }
  )

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    if (requireNamespace("cli", quietly = TRUE)) cli::cli_alert_info(paste("Testing for", k_cuts, "cut-point(s)..."))

    res <- .systematic_search(userdata, k_cuts, "p_value", covariates, nmin, "factor", quiet = TRUE)
    if (any(is.na(res$optimal_cuts))) {
      if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("No valid cut-points found for", k_cuts, "cut(s) due to model failures or constraints."))
      results <- rbind(results, data.frame(num_cuts = k_cuts, IC = NA_real_, cuts = I(list(NULL))))
      next
    }
    factor_status <- factor(findInterval(userdata$factor, res$optimal_cuts, left.open = TRUE) + 1L)
    best_ic <- .get_model_ic_num(userdata, factor_status, k_cuts, n, criterion, cov_part)

    results <- rbind(results, data.frame(num_cuts = k_cuts, IC = best_ic, cuts = I(list(res$optimal_cuts))))
  }
  names(results)[2] <- criterion
  return(results)
}

#' Internal helper: Compute IC from standard factors
#'
#' @description
#' Computes the Information Criterion dynamically for the systematic engine
#' based on the best discovered grouping factors.
#'
#' @param userdata Cleaned survival data frame.
#' @param factor_status The cut survival factor.
#' @param k_cuts Number of cuts evaluated.
#' @param n Sample size.
#' @param criterion IC to calculate (AIC, AICc, BIC).
#' @param cov_part Formula string for covariates.
#'
#' @return Single numeric IC value.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Computes AIC, AICc, or BIC from a fitted Cox model.
#' @srrstats {RE1.0} Relies on `logLik()` from `coxph` fit.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @noRd
.get_model_ic_num <- function(userdata, factor_status, k_cuts, n, criterion, cov_part) {
  num_cov <- length(cov_part[cov_part != ""])
  formula_str <- paste("survival::Surv(time, event) ~ factor_status", cov_part)
  fit <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = userdata), error = function(e) NULL)
  if (is.null(fit)) {
    return(Inf)
  }
  .calc_ic(fit, k_cuts + num_cov, n, criterion)
}
