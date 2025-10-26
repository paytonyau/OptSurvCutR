#' Validate an Optimal Cut-point Using Bootstrapping
#'
#' @description
#' Assesses the stability of optimal cut-points found by \code{\link{find_cutpoint}} by
#' performing a bootstrap analysis and generating 95\% confidence intervals.
#' This version is streamlined for survival (time-to-event) analysis.
#'
#' @param cutpoint_result An object returned by the \code{\link{find_cutpoint}} function.
#' @param num_replicates The number of bootstrap replicates to perform. Default is 500.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores for bootstrapping.
#' @param n_cores Number of cores to use for parallel processing; if NULL, defaults to 2.
#' @param seed An optional integer to set the random seed for reproducible results.
#' @param nmin The minimum group size to enforce during bootstrap runs. Defaults
#'   to 90\% of the `nmin` from the original analysis to reduce failures
#'   in sparse bootstrap samples.
#' @param ... Additional arguments to be passed to \code{\link{find_cutpoint}}, especially
#'   for the genetic algorithm (e.g., `popSize`, `maxiter`).
#'
#' @return An object of class `validate_cutpoint_result` containing the
#'   original cut-point(s), their 95\% confidence intervals, bootstrap distribution,
#'   and analysis parameters.
#' @importFrom foreach %dopar%
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_alert_warning cli_inform cli_abort cli_progress_bar cli_progress_update
#' @importFrom stats quantile na.omit complete.cases sd median IQR
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom tidyr pivot_longer everything
#' @export
validate_cutpoint <- function(cutpoint_result, num_replicates = 500, use_parallel = FALSE, n_cores = NULL, seed = NULL, nmin = NULL, ...) {
  # --- 1. Validate Input and Set Seed ---
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    cli::cli_abort("Input must be an object from the find_cutpoint function.")
  }
  if (!is.numeric(num_replicates) || num_replicates < 1 || num_replicates != round(num_replicates)) {
    cli::cli_abort("num_replicates must be a positive integer.")
  }
  if (num_replicates < 20) {
    cli::cli_abort("num_replicates must be at least 20 for reliable validation.")
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

  if (any(is.na(original_cuts))) {
    cli::cli_abort("Input 'find_cutpoint' object contains NA cut-points. Cannot validate.")
  }

  # Extract predictor and outcome variables
  predictor <- original_params$predictor
  num_cuts <- original_params$num_cuts
  method <- original_params$method
  criterion <- original_params$criterion
  covariates <- original_params$covariates

  # Sensible default for nmin in bootstrap samples
  if (is.null(nmin)) {
    nmin <- floor(0.9 * original_params$nmin)
    cli::cli_alert_info("`nmin` for bootstrap not specified, defaulting to 90% of original: {nmin}")
  }
  if (!is.numeric(nmin) || nmin <= 0) {
    cli::cli_abort("nmin must be a positive number.")
  }
  if (n < nmin * (num_cuts + 1)) {
    cli::cli_abort("Not enough data ({n}) for nmin ({nmin}) and {num_cuts} cut(s).")
  }

  cli::cli_alert_info("Validating {num_cuts} cut-point(s) from '{method}' search using '{criterion}' criterion.")

  # --- 3. Setup Parallel Backend ---
  cores_to_use <- if (use_parallel) {
    cores_available <- parallel::detectCores() %||% 2
    min(cores_available - 1, n_cores %||% 2, num_replicates)
  } else {
    1
  }
  if (cores_to_use < 1) cores_to_use <- 1

  if (use_parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      cli::cli_abort("Package 'doParallel' is required for parallel processing.")
    }
    cl <- parallel::makeCluster(cores_to_use)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cli::cli_alert_info("Running {num_replicates} bootstrap replicates on {cores_to_use} core{?s}...")
  } else {
    foreach::registerDoSEQ()
    cli::cli_alert_info("Running {num_replicates} bootstrap replicates sequentially...")
  }

  # --- 4. Main Bootstrap Loop ---
  # Initialize progress bar for sequential execution
  if (!use_parallel) {
    pb <- cli::cli_progress_bar("Bootstrapping", total = num_replicates)
  }

  # Export necessary functions for parallel processing
  functions_to_export <- if (method == "genetic") {
    c("find_cutpoint", ".run_genetic_search", ".systematic_search", ".get_stat")
  } else {
    c("find_cutpoint", ".systematic_search", ".get_stat")
  }

  i <- NULL

  bootstrap_results <- foreach::foreach(
    i = 1:num_replicates,
    .combine = 'rbind',
    .export = functions_to_export,
    .errorhandling = 'pass'
  ) %dopar% {
    # Update progress bar for sequential runs
    if (!use_parallel) {
      cli::cli_progress_update(id = pb)
    }

    # Set seed per replicate for reproducibility
    set.seed(i)

    # Resample with replacement
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- original_data[boot_indices, ]

    # Prepare arguments for find_cutpoint
    extra_args <- list(...)
    final_args <- c(
      list(
        data = boot_data,
        predictor = "factor",  # Matches internal name in find_cutpoint
        outcome_time = "time",
        outcome_event = "event",
        method = method,
        num_cuts = num_cuts,
        criterion = criterion,
        covariates = covariates,
        nmin = nmin,
        quiet = TRUE
      ),
      extra_args
    )

    # Run find_cutpoint on bootstrap sample
    res <- tryCatch({
      suppressMessages(do.call(find_cutpoint, final_args))
    }, error = function(e) {
      cli::cli_inform("Bootstrap replicate {i} failed: {e$message}")
      return(NULL)
    })

    if (is.null(res) || any(is.na(res$optimal_cuts)) || length(res$optimal_cuts) != num_cuts) {
      return(rep(NA, num_cuts))
    }
    res$optimal_cuts
  }

  # --- 5. Process Results and Calculate CIs ---
  # Ensure bootstrap_results is a matrix
  if (!is.matrix(bootstrap_results)) {
    bootstrap_matrix <- matrix(bootstrap_results, nrow = num_replicates, ncol = num_cuts, byrow = TRUE)
  } else {
    bootstrap_matrix <- bootstrap_results
  }
  colnames(bootstrap_matrix) <- paste0("Cut", 1:num_cuts)

  # Check for failed replicates
  successful_reps <- sum(stats::complete.cases(bootstrap_matrix))
  failed_reps <- num_replicates - successful_reps

  if (failed_reps > 0) {
    cli::cli_alert_warning("{failed_reps} of {num_replicates} bootstrap replicates failed to find valid cut-points.")
  }
  if (successful_reps < 20) {
    cli::cli_abort("Fewer than 20 successful bootstrap replicates ({successful_reps}) completed. Increase num_replicates or reduce nmin.")
  }

  cli::cli_alert_success("{successful_reps} replicates completed successfully.")

  # Warning for low success rate
  warn_threshold <- 0.8
  if (successful_reps / num_replicates < warn_threshold && successful_reps > 0) {
    warning(
      paste(successful_reps, "of", num_replicates, "bootstrap replicates succeeded.",
            "The confidence intervals may be unreliable.")
    )
  }

  # Clean bootstrap results and calculate statistics
  bootstrap_matrix_clean <- na.omit(bootstrap_matrix)

  # Calculate 95% CI for each cut-point
  ci <- apply(bootstrap_matrix_clean, 2, stats::quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  if (is.vector(ci)) {
    ci_df <- data.frame(Lower = ci[1], Upper = ci[2])
  } else {
    ci_df <- as.data.frame(t(ci))
    names(ci_df) <- c("Lower", "Upper")
  }
  row.names(ci_df) <- paste0("Cut ", 1:num_cuts)

  # Calculate descriptive statistics
  boot_summary <- lapply(1:num_cuts, function(i) {
    valid_cuts <- bootstrap_matrix_clean[, i]
    if (length(valid_cuts) == 0) {
      return(list(mean = NA, sd = NA, median = NA, Q1 = NA, Q3 = NA))
    }
    list(
      mean = mean(valid_cuts, na.rm = TRUE),
      sd = stats::sd(valid_cuts, na.rm = TRUE),
      median = stats::median(valid_cuts, na.rm = TRUE),
      Q1 = stats::quantile(valid_cuts, 0.25, na.rm = TRUE),
      Q3 = stats::quantile(valid_cuts, 0.75, na.rm = TRUE)
    )
  })
  names(boot_summary) <- paste0("Cut", 1:num_cuts)

  # Prepare bootstrap distribution as data frame
  bootstrap_df <- as.data.frame(bootstrap_matrix_clean)
  names(bootstrap_df) <- paste0("Cut_point_", 1:num_cuts)

  # Construct output object
  output <- list(
    original_cuts = original_cuts,
    confidence_intervals = ci_df,
    bootstrap_distribution = bootstrap_df,
    boot_summary = boot_summary,
    parameters = list(
      num_replicates = num_replicates,
      successful_reps = successful_reps,
      failed_reps = failed_reps,
      use_parallel = use_parallel,
      n_cores = cores_to_use,
      seed = seed,
      nmin = nmin,
      method = method,
      criterion = criterion,
      covariates = covariates
    )
  )
  class(output) <- "validate_cutpoint_result"

  print(output)
  invisible(output)
}

#' @param x An object of class `validate_cutpoint_result`.
#' @param ... Unused.
#' @rdname validate_cutpoint
#' @export
print.validate_cutpoint_result <- function(x, ...) {
  cat("Cut-point Stability Analysis (Bootstrap)\n")
  cat("----------------------------------------\n")
  cat("Original Optimal Cut-point(s):", paste(round(x$original_cuts, 3), collapse = ", "), "\n")
  cat("Successful Replicates:", x$parameters$successful_reps, "/", x$parameters$num_replicates,
      "(", round(100 * x$parameters$successful_reps / x$parameters$num_replicates, 1), "%)\n")
  cat("Failed Replicates:", x$parameters$failed_reps, "\n\n")

  cat("95% Confidence Intervals\n")
  cat("------------------------\n")
  print(round(x$confidence_intervals, 3))

  cat("\nBootstrap Summary Statistics\n")
  cat("---------------------------\n")
  summary_df <- do.call(rbind, lapply(names(x$boot_summary), function(cut) {
    stats <- x$boot_summary[[cut]]
    data.frame(
      Cut = cut,
      Mean = stats$mean,
      SD = stats$sd,
      Median = stats$median,
      Q1 = stats$Q1,
      Q3 = stats$Q3
    )
  }))
  numeric_cols <- sapply(summary_df, is.numeric)
  summary_df[, numeric_cols] <- round(summary_df[, numeric_cols], 3)
  print(summary_df)

  cat("\nHint: Use `summary()` for detailed statistics or `plot()` to visualize stability.\n")
  invisible(x)
}

#' @rdname validate_cutpoint
#' @importFrom ggplot2 ggplot aes .data geom_density geom_vline labs theme_minimal facet_wrap
#' @export
plot.validate_cutpoint_result <- function(x, ...) {
  dist_data <- x$bootstrap_distribution
  num_cuts <- ncol(dist_data)

  if (x$parameters$successful_reps == 0) {
    cli::cli_inform("Cannot generate plot due to no successful bootstrap replicates.")
    return(invisible(NULL))
  }

  plot_data <- tidyr::pivot_longer(
    dist_data,
    cols = tidyr::everything(),
    names_to = "Cut",
    values_to = "Value",
    names_prefix = "Cut_point_"
  )
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
    ggplot2::labs(
      title = "Bootstrap Distribution of Optimal Cut-points",
      subtitle = paste(x$parameters$successful_reps, "successful replicates"),
      x = "Cut-point Value",
      y = "Density"
    ) +
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
  cat("Cut-point Stability Analysis (Bootstrap)\n")
  cat("----------------------------------------\n")
  cat("Original Optimal Cut-point(s):", paste(round(object$original_cuts, 3), collapse = ", "), "\n\n")

  if (show_descriptives) {
    cat("Bootstrap Distribution Summary\n")
    cat("-----------------------------\n")
    summary_df <- do.call(rbind, lapply(names(object$boot_summary), function(cut) {
      stats <- object$boot_summary[[cut]]
      data.frame(
        Cut = cut,
        Mean = stats$mean,
        SD = stats$sd,
        Median = stats$median,
        Q1 = stats$Q1,
        Q3 = stats$Q3
      )
    }))
    numeric_cols <- sapply(summary_df, is.numeric)
    summary_df[, numeric_cols] <- round(summary_df[, numeric_cols], 3)
    print(summary_df)
    cat("\n")
  }

  if (show_ci) {
    cat("95% Confidence Intervals\n")
    cat("------------------------\n")
    print(round(object$confidence_intervals, 3))
    cat("\n")
  }

  if (show_params) {
    cat("Validation Parameters\n")
    cat("---------------------\n")
    cat("Replicates Requested:", object$parameters$num_replicates, "\n")
    cat("Successful Replicates:", object$parameters$successful_reps,
        "(", round(100 * object$parameters$successful_reps / object$parameters$num_replicates, 1), "%)\n")
    cat("Failed Replicates:", object$parameters$failed_reps, "\n")
    cat("Parallel Processing:", ifelse(object$parameters$use_parallel, "Enabled", "Disabled"), "\n")
    cat("Cores Used:", object$parameters$n_cores, "\n")
    cat("Seed:", ifelse(is.null(object$parameters$seed), "Not set", object$parameters$seed), "\n")
    cat("Minimum Group Size (nmin):", object$parameters$nmin, "\n")
    cat("Method:", object$parameters$method, "\n")
    cat("Criterion:", object$parameters$criterion, "\n")
    cat("Covariates:", ifelse(is.null(object$parameters$covariates), "None", paste(object$parameters$covariates, collapse = ", ")), "\n")
    cat("\n")
  }

  if (plot.it) {
    cat("Bootstrap Distribution Plot\n")
    cat("--------------------------\n")
    print(plot(object, ...))
  }

  invisible(object)
}
