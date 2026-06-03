# ===================================================================
# INTERNAL UTILITY: PERMUTATION TESTING
# Calculates exact p-values to correct for optimization bias
# ===================================================================

#' Internal helper: Permutation Testing for Optimization Bias
#'
#' @description
#' Performs permutation testing by shuffling survival outcomes (time and event)
#' as a locked pair to construct a null distribution of the maximally selected rank statistic.
#' This calculates an empirically adjusted p-value, correcting for the Type I
#' error inflation inherently caused by multi-dimensional threshold searching.
#'
#' @param time_vec Numeric vector. Survival time values.
#' @param censor_vec Numeric vector. Survival event status values (0 or 1).
#' @param predictor_vec Numeric vector. Continuous predictor biomarker values.
#' @param num_cuts Integer. Number of cut-points being evaluated.
#' @param method Character. Search algorithm ("systematic" or "genetic").
#' @param criterion Character. Statistic being optimized.
#' @param covariates Character vector. Optional covariates.
#' @param userdata_full Cleaned data frame template containing covariate metrics.
#' @param nmin Integer. Minimum observations per group.
#' @param obs_stat Numeric. The observed maximum statistic from the true data.
#' @param n_perm Integer. Number of permutations to run.
#' @param n_cores Integer. Number of CPU cores for parallel processing.
#' @param max.generations Integer. Maximum generations for genetic search.
#' @param pop.size Integer. Population size for genetic search.
#' @param use_cpp Logical. If TRUE, deploys the optimized Rcpp matrix factory backend.
#' @param grid_by Numeric. Percentile step increment for systematic downsampling.
#' @param ... Additional arguments passed to the core engines.
#'
#' @return A single numeric value representing the permuted p-value, or `NA_real_`
#' if the calculation fails.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.3} Computes empirical p-values via permutation testing to account for multiple comparisons.
#' @srrstats {G5.2} Gracefully handles failed permutation replicates by dropping them before calculation.
#' @srrstats {G2.14a} Handles `NA`s generated during resampling loops using `na.omit()`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom foreach %dopar% foreach registerDoSEQ
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel makeCluster stopCluster
#' @importFrom stats na.omit
#' @noRd
.run_permutations <- function(time_vec, censor_vec, predictor_vec, num_cuts,
                              method, criterion, covariates, userdata_full, nmin,
                              obs_stat, n_perm, n_cores, max.generations, pop.size,
                              use_cpp = TRUE, grid_by = 0.01, ...) {

  # Set up parallel or sequential backend
  if (n_cores > 1 && requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    # Ensure cluster closes even if function crashes (prevents zombie processes)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  n_obs <- length(time_vec)
  extra_args <- list(...)

  i <- NULL
  null_stats <- foreach::foreach(
    i = 1:n_perm,
    .combine = c,
    .errorhandling = "remove", # Silently drops runs that fail due to extreme shuffles
    .packages = c("survival", "OptSurvCutR")
  ) %dopar% {

    # SAFE SHUFFLING TASK: Shuffle outcomes as an intact locked pair to protect censoring histories
    shuffle_idx <- sample(n_obs)
    shuffled_time  <- time_vec[shuffle_idx]
    shuffled_event <- censor_vec[shuffle_idx]

    # Reconstruct the tracking dataset iteration shell cleanly without altering baseline covariates
    shuffled_data <- userdata_full
    shuffled_data$time <- shuffled_time
    shuffled_data$event <- shuffled_event

    # Route to the appropriate accelerated engine
    if (method == "systematic") {
      res <- do.call(.systematic_search, c(
        list(userdata = shuffled_data, num_cuts = num_cuts, criterion = criterion,
             covariates = covariates, nmin = nmin, predictor_name = "factor",
             use_cpp = use_cpp, grid_by = grid_by, quiet = TRUE),
        extra_args
      ))
      return(res$optimal_stat)
    } else {
      confound_df <- if (!is.null(covariates)) shuffled_data[, covariates, drop = FALSE] else NULL
      res <- do.call(.run_genetic_search, c(
        list(target = predictor_vec, numcut = num_cuts, time = shuffled_time,
             censor = shuffled_event, confound = confound_df, nmin = nmin,
             criterion = criterion, max.generations = max.generations,
             pop.size = pop.size, use_cpp = use_cpp, print.level = 0),
        extra_args
      ))
      if (is.null(res) || !is.finite(res$value)) return(NA_real_)
      return(res$value)
    }
  }

  # Calculate empirical p-value
  valid_nulls <- stats::na.omit(null_stats)
  if (length(valid_nulls) == 0) return(NA_real_)

  p_perm <- (sum(valid_nulls >= obs_stat) + 1) / (length(valid_nulls) + 1)
  return(p_perm)
}
