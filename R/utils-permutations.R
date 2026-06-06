# ===================================================================
# INTERNAL UTILITY: PERMUTATION TESTING
# Calculates exact p-values to correct for optimization bias
# ===================================================================

#' Internal helper: Permutation Testing for Optimization Bias
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
  n_obs <- length(time_vec)
  extra_args <- list(...)
  i <- NULL

  # Dynamically configure error handling: "pass" keeps covr tracking alive,
  # while "remove" stays active for standard production performance.
  is_covr <- Sys.getenv("R_COVR") == "true"
  err_handle <- if (is_covr) "pass" else "remove"

  # Enforce foreground sequential execution during code coverage runs
  if (n_cores > 1 && !is_covr && requireNamespace("doParallel", quietly = TRUE)) {
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
  } else {
    foreach::registerDoSEQ()
  }

  null_stats <- foreach::foreach(
    i = 1:n_perm,
    .combine = c,
    .errorhandling = err_handle,
    .packages = c("survival", "OptSurvCutR")
  ) %dopar% {
    shuffle_idx <- sample(n_obs)
    shuffled_time <- time_vec[shuffle_idx]
    shuffled_event <- censor_vec[shuffle_idx]

    shuffled_data <- userdata_full
    shuffled_data$time <- shuffled_time
    shuffled_data$event <- shuffled_event

    if (method == "systematic") {
      res <- do.call(.systematic_search, c(
        list(
          userdata = shuffled_data, num_cuts = num_cuts, criterion = criterion,
          covariates = covariates, nmin = nmin, predictor_name = "factor",
          use_cpp = use_cpp, grid_by = grid_by, quiet = TRUE
        ),
        extra_args
      ))
      return(res$optimal_stat)
    } else {
      confound_df <- if (!is.null(covariates)) shuffled_data[, covariates, drop = FALSE] else NULL

      res <- do.call(.run_genetic_search, c(
        list(
          target = predictor_vec, numcut = num_cuts, time = shuffled_time,
          censor = shuffled_event, confound = confound_df, nmin = nmin,
          criterion = criterion, max.generations = max.generations,
          pop.size = pop.size, use_cpp = use_cpp, print.level = 0
        ),
        extra_args
      ))

      if (is.null(res) || !is.finite(res$value)) {
        return(NA_real_)
      }
      return(res$value)
    }
  }

  # Cleanly filter out error structures or NA values generated during optimization drops
  valid_nulls <- na.omit(sapply(null_stats, function(x) {
    if (is.numeric(x) && is.finite(x)) {
      return(x)
    } else {
      return(NA_real_)
    }
  }))

  if (length(valid_nulls) == 0) {
    return(NA_real_)
  }

  p_perm <- (sum(valid_nulls >= obs_stat) + 1) / (length(valid_nulls) + 1)
  return(p_perm)
}
