#' Validate an Optimal Cut-point Using Bootstrapping
#'
#' @description
#' Assesses the stability of optimal cut-points found by `find_cutpoint` by
#' performing a bootstrap analysis and generating 95% confidence intervals.
#' This version is streamlined for survival (time-to-event) analysis.
#'
#' @param cutpoint_result An object returned by the `find_cutpoint` function.
#' @param num_replicates The number of bootstrap replicates to perform. Default is 500.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores for bootstrapping.
#' @param seed An optional integer to set the random seed for reproducible results.
#' @param nmin The minimum group size to enforce during bootstrap runs. Defaults
#'   to 90% of the `nmin` from the original analysis to reduce failures
#'   in sparse bootstrap samples.
#' @param ... Additional arguments to be passed to `find_cutpoint`, especially
#'   for the genetic algorithm (e.g., `popSize`, `maxiter`).
#'
#' @return An object of class `validate_cutpoint_result` containing the
#'   original cut-point(s) and their 95% confidence intervals.
#' @importFrom foreach %dopar%
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_bullets cli_warn cli_progress_bar cli_progress_update
#' @importFrom stats quantile na.omit complete.cases sd median IQR
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom tidyr pivot_longer everything
#' @export
validate_cutpoint <- function(cutpoint_result, num_replicates = 500, use_parallel = FALSE, seed = NULL, nmin = NULL, ...) {

  # --- 1. Validate Input and Set Seed ---
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    stop("Input must be an object from the find_cutpoint function.", call. = FALSE)
  }
  if (!is.null(seed)) {
    set.seed(seed)
    cli::cli_alert_info("Using random seed: {seed} for reproducible bootstrapping.")
  }

  # --- 2. Extract Original Parameters and Data ---
  original_params <- cutpoint_result$parameters
  original_data <- cutpoint_result$userdata
  original_cuts <- cutpoint_result$optimal_cuts
  n <- nrow(original_data)

  # Sensible default for nmin in bootstrap samples
  if (is.null(nmin)) {
    nmin <- floor(0.9 * original_params$nmin)
    cli::cli_alert_info("`nmin` for bootstrap not specified, defaulting to 90% of original: {nmin}")
  }

  cli::cli_alert_info("Validating {length(original_cuts)} cut-point(s) from '{original_params$method}' search.")

  # --- 3. Setup Parallel Backend ---
  if (use_parallel) {
    cores_available <- parallel::detectCores()
    cores_to_use <- min(cores_available, num_replicates)
    cl <- parallel::makeCluster(cores_to_use)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cli::cli_alert_info("Running {num_replicates} bootstrap replicates on {cores_to_use} cores...")
  } else {
    foreach::registerDoSEQ()
    cli::cli_alert_info("Running {num_replicates} bootstrap replicates sequentially...")
  }

  # --- 4. Main Bootstrap Loop ---

  # Initialize progress bar for sequential execution
  if (!use_parallel) {
    pb <- cli::cli_progress_bar("Bootstrapping", total = num_replicates)
  }

  # We no longer need to manually export `%||%` because it's an exported function
  functions_to_export <- c("find_cutpoint", ".systematic_search", ".get_stat",
                           ".run_genetic_search", ".obj")

  bootstrap_results <- foreach::foreach(
    i = 1:num_replicates,
    .combine = 'rbind',
    .export = functions_to_export
  ) %dopar% {

    # Update progress bar inside the loop ONLY for sequential runs
    if (!use_parallel) {
      cli::cli_progress_update(id = pb)
    }

    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- original_data[boot_indices, ]

    extra_args <- list(...)

    final_args <- c(
      list(
        data = boot_data,
        predictor = "factor",
        outcome_time = "time",
        outcome_event = "event",
        method = original_params$method,
        num_cuts = original_params$num_cuts,
        criterion = original_params$criterion,
        covariates = original_params$covariates,
        nmin = nmin,
        quiet = TRUE
      ),
      extra_args
    )

    res <- tryCatch({
      suppressMessages(do.call(find_cutpoint, final_args))
    }, error = function(e) NULL)

    if (is.null(res) || any(is.na(res$optimal_cuts)) || length(res$optimal_cuts) != length(original_cuts)) {
      return(rep(NA, length(original_cuts)))
    } else {
      return(res$optimal_cuts)
    }
  }

  # --- 5. Process Results and Calculate CIs ---
  if (!is.matrix(bootstrap_results)) {
    bootstrap_matrix <- as.matrix(t(bootstrap_results))
  } else {
    bootstrap_matrix <- bootstrap_results
  }

  successful_reps <- sum(stats::complete.cases(bootstrap_matrix))
  failed_reps <- num_replicates - successful_reps

  if (failed_reps > 0) {
    cli::cli_warn("{failed_reps} of {num_replicates} bootstrap replicates failed to find a valid cut-point.")
  }
  cli::cli_alert_success("{successful_reps} replicates completed successfully.")

  min_required_reps <- 20
  if (successful_reps < min_required_reps) {
    stop(paste0("Fewer than ", min_required_reps, " successful bootstrap replicates completed. ",
                "Cannot calculate a reliable confidence interval."), call. = FALSE)
  }

  bootstrap_matrix_clean <- na.omit(bootstrap_matrix)

  ci <- apply(bootstrap_matrix_clean, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)

  if (is.vector(ci)) {
    ci_df <- data.frame(Lower = ci[1], Upper = ci[2])
  } else {
    ci_df <- as.data.frame(t(ci))
    names(ci_df) <- c("Lower", "Upper")
  }
  row.names(ci_df) <- paste0("Cut ", 1:length(original_cuts))

  bootstrap_df <- as.data.frame(bootstrap_matrix_clean)
  names(bootstrap_df) <- paste0("Cut_point_", 1:ncol(bootstrap_df))

  output <- list(original_cuts = original_cuts,
                 confidence_intervals = ci_df,
                 bootstrap_distribution = bootstrap_df,
                 parameters = list(num_replicates = num_replicates,
                                   successful_reps = successful_reps))
  class(output) <- "validate_cutpoint_result"

  print(output)
  invisible(output)
}


#' @param x An object of class `validate_cutpoint_result`.
#' @param ... Unused.
#' @rdname validate_cutpoint
#' @export
print.validate_cutpoint_result <- function(x, ...) {
  cli::cli_h1("Cut-point Stability Analysis (Bootstrap)")
  cli::cli_bullets(c(
    "*" = "Original Optimal Cut-point(s): {.strong {paste(round(x$original_cuts, 3), collapse=', ')}}",
    "*" = "Successful Replicates: {x$parameters$successful_reps} / {x$parameters$num_replicates}"
  ))
  cli::cli_text("\n95% Confidence Intervals:")
  print(round(x$confidence_intervals, 3))
  cli::cli_alert_info("\nHint: Use `summary()` for detailed statistics or `plot()` to visualize stability.")
}


#' @rdname validate_cutpoint
#' @importFrom ggplot2 ggplot aes .data geom_density geom_vline labs theme_minimal facet_wrap
#' @export
plot.validate_cutpoint_result <- function(x, ...) {

  dist_data <- x$bootstrap_distribution
  num_cuts <- ncol(dist_data)

  plot_data <- tidyr::pivot_longer(dist_data,
                                   cols = tidyr::everything(),
                                   names_to = "Cut",
                                   values_to = "Value",
                                   names_prefix = "Cut_point_")
  plot_data$Cut <- paste("Cut-point", plot_data$Cut)

  line_data <- data.frame(
    Cut = paste("Cut-point", 1:num_cuts),
    original_cut = x$original_cuts,
    ci_lower = x$confidence_intervals$Lower,
    ci_upper = x$confidence_intervals$Upper
  )

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Value)) +
    ggplot2::geom_density(fill = "#56B4E9", color = "#0072B2", alpha = 0.6) +
    ggplot2::geom_vline(
      data = line_data,
      aes(xintercept = .data$original_cut),
      color = "#D55E00", linetype = "solid", linewidth = 1
    ) +
    ggplot2::geom_vline(
      data = line_data,
      aes(xintercept = .data$ci_lower),
      color = "#D55E00", linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      data = line_data,
      aes(xintercept = .data$ci_upper),
      color = "#D55E00", linetype = "dashed"
    ) +
    ggplot2::labs(title = "Bootstrap Distribution of Optimal Cut-points",
                  subtitle = paste(x$parameters$successful_reps, "successful replicates"),
                  x = "Cut-point Value", y = "Density") +
    ggplot2::theme_minimal(base_size = 14)

  if (num_cuts > 1) {
    p <- p + ggplot2::facet_wrap(~Cut, scales = "free_x")
  }

  return(p)
}

#' @param object An object of class `validate_cutpoint_result`.
#' @param show_descriptives Logical. If TRUE, shows descriptive statistics of the bootstrap distributions.
#' @param show_ci Logical. If TRUE, shows the confidence intervals.
#' @param show_params Logical. If TRUE, shows the parameters of the validation run.
#' @param plot.it Logical. If TRUE, displays the density plot of the bootstrap distributions.
#' @rdname validate_cutpoint
#' @export
summary.validate_cutpoint_result <- function(object, show_descriptives = TRUE, show_ci = TRUE, show_params = TRUE, plot.it = FALSE, ...) {

  cli::cli_h1("Cut-point Stability Analysis (Bootstrap)")
  cli::cli_text("Original Optimal Cut-point(s): {.strong {paste(round(object$original_cuts, 3), collapse=', ')}}")

  if (show_descriptives) {
    cli::cli_h2("Bootstrap Distribution Summary")

    dist_data <- object$bootstrap_distribution

    desc_stats <- data.frame(
      Mean = sapply(dist_data, mean, na.rm = TRUE),
      Median = sapply(dist_data, stats::median, na.rm = TRUE),
      SD = sapply(dist_data, stats::sd, na.rm = TRUE),
      IQR = sapply(dist_data, stats::IQR, na.rm = TRUE)
    )
    row.names(desc_stats) <- paste("Cut", 1:nrow(desc_stats))

    print(round(desc_stats, 3))
  }

  if (show_ci) {
    cli::cli_h2("Confidence Intervals")
    print(round(object$confidence_intervals, 3))
  }

  if (show_params) {
    cli::cli_h2("Validation Parameters")
    cli::cli_bullets(c(
      "*" = "Replicates Requested: {object$parameters$num_replicates}",
      "*" = "Successful Replicates: {object$parameters$successful_reps} ({round(100 * object$parameters$successful_reps / object$parameters$num_replicates, 1)}%)"
    ))
  }

  if (plot.it) {
    cli::cli_h2("Bootstrap Distribution Plot")
    print(plot(object, ...))
  }

  invisible(object)
}


# Suppress NOTE about ggplot2/tidyr global variables
utils::globalVariables(c("Value", "everything"))
