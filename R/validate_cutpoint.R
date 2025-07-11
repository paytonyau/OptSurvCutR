#' Validate an Optimal Cut-point Using Bootstrapping
#'
#' Assesses the stability of optimal cut-points found by `find_cutpoint` by
#' performing a bootstrap analysis and generating 95% confidence intervals.
#'
#' @param cutpoint_result An object returned by the `find_cutpoint` function.
#' @param num_replicates The number of bootstrap replicates to perform.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores for bootstrapping.
#' @param nmin The minimum group size to use during bootstrap runs.
#'
#' @return An object of class `validate_cutpoint_result` containing the
#'   original cut-point(s) and their 95% confidence intervals.
#' @importFrom foreach %dopar%
#' @export
validate_cutpoint <- function(cutpoint_result, num_replicates = 500, use_parallel = FALSE, nmin = 15) {

  # --- 1. Validate Input and Extract Original Data ---
  if (!inherits(cutpoint_result, c("find_cutpoint_systematic", "find_cutpoint_genetic"))) {
    stop("Input must be an object from the find_cutpoint function.", call. = FALSE)
  }

  original_data <- cutpoint_result$userdata
  original_params <- cutpoint_result$parameters
  n <- nrow(original_data)

  if (inherits(cutpoint_result, "find_cutpoint_systematic")) {
    original_cuts <- cutpoint_result$best_by_vote
    cli::cli_alert_info("Validating 1 cut-point from systematic search.")
  } else {
    original_cuts <- cutpoint_result$optimal_cuts
    cli::cli_alert_info("Validating {length(original_cuts)} cut-point(s) from genetic algorithm.")
  }

  # --- 2. Setup Parallel Backend ---
  if (use_parallel) {
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cli::cli_alert_info("Running {num_replicates} bootstrap replicates on {cores} cores...")
  } else {
    foreach::registerDoSEQ()
    cli::cli_alert_info("Running {num_replicates} bootstrap replicates sequentially...")
  }

  # --- 3. Main Bootstrap Loop ---
  bootstrap_results <- foreach::foreach(i = 1:num_replicates, .combine = 'rbind') %dopar% {
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- original_data[boot_indices, ]

    original_method <- if(inherits(cutpoint_result, "find_cutpoint_systematic")) "systematic" else "genetic"

    base_args <- list(
      data = boot_data,
      predictor = "factor",
      method = original_method,
      num_cuts = length(original_cuts),
      nmin = nmin,
      use_parallel = FALSE
    )

    if (original_params$analysis_type == "survival") {
      outcome_args <- list(outcome_time = "time", outcome_event = "event")
    } else {
      outcome_args <- list(outcome_binary = "outcome")
    }

    final_args <- c(base_args, outcome_args)

    res <- tryCatch({
      do.call(find_cutpoint, final_args)
    }, error = function(e) NULL)

    if (is.null(res)) { return(NA)
    } else {
      if (inherits(res, "find_cutpoint_systematic")) { return(res$best_by_vote)
      } else { return(res$optimal_cuts) }
    }
  }

  bootstrap_matrix <- na.omit(bootstrap_results)
  successful_reps <- nrow(bootstrap_matrix)
  cli::cli_alert_success("{successful_reps} of {num_replicates} bootstrap replicates completed successfully.")

  # --- START OF CORRECTED SECTION ---
  # Check if at least 80% of replicates succeeded, with a minimum of 10.
  min_required_reps <- max(10, floor(0.8 * num_replicates))
  if (successful_reps < min_required_reps) {
    stop(paste0("Too few successful bootstrap replicates (", successful_reps,
                ") to calculate a reliable confidence interval. Try using a smaller 'nmin'."),
         call. = FALSE)
  }
  # --- END OF CORRECTED SECTION ---

  ci <- apply(bootstrap_matrix, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  if (is.vector(ci)) {
    ci_df <- data.frame(Lower = ci[1], Upper = ci[2])
  } else {
    ci_df <- as.data.frame(t(ci))
    names(ci_df) <- c("Lower", "Upper")
  }
  row.names(ci_df) <- paste0("Cut ", 1:length(original_cuts))

  output <- list(original_cuts = original_cuts, confidence_intervals = ci_df,
                 bootstrap_distribution = bootstrap_matrix)
  class(output) <- "validate_cutpoint_result"

  print(output)
  invisible(output)
}


#' Print method for validate_cutpoint_result
#' @param x An object of class `validate_cutpoint_result`.
#' @param ... Unused.
#' @export
print.validate_cutpoint_result <- function(x, ...) {
  cli::cli_h1("Cut-point Stability Analysis (Bootstrap)")
  cli::cli_bullets(c("Original Optimal Cut-point(s):" = paste(round(x$original_cuts, 3), collapse=", ")))
  cli::cli_text("\n95% Confidence Intervals:")
  print(round(x$confidence_intervals, 3))
  cli::cli_alert_info("\nHint: Narrow confidence intervals suggest stable and reliable cut-points.")
}
