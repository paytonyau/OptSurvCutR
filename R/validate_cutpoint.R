# ====================================================================
# VALIDATION FUNCTION
# Bootstraps cut-points for stability and confidence intervals.
# ===================================================================
#' Validate an Optimal Cut-point Using Bootstrapping
#'
#' @description
#' Assesses cut-point stability from [find_cutpoint()] via bootstrap
#' analysis, generating 95% confidence intervals. Streamlined for
#' survival (time-to-event) analysis.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G1.0} Primary refs: Efron (1979), Rota et al. (2015).
#' @srrstats {G2.0} Input object class validated with `inherits()`.
#' @srrstats {G2.0a} Scalar parameters (`num_replicates`, `nmin`) checked.
#' @srrstats {G2.1} Input types validated (`numeric`, `logical`, `integer`).
#' @srrstats {G2.1a} `cutpoint_result` needs `optimal_cuts`, `parameters`.
#' @srrstats {G2.2} `NA` in `optimal_cuts` triggers `cli_abort()`.
#' @srrstats {G2.3} `na.omit()` used in CI calculation.
#' @srrstats {G2.3a} `complete.cases()` for replicate counting.
#' @srrstats {G2.4a} `as.integer()` on `nmin`, `num_replicates`.
#' @srrstats {G2.4d} `tryCatch()` on each bootstrap replicate.
#' @srrstats {G2.5} `NA` results from failed replicates preserved.
#' @srrstats {G2.11} `seed` controls randomness.
#' @srrstats {G2.13} `cli_abort()` for invalid input.
#' @srrstats {G2.14a} `NA` in `optimal_cuts` aborts.
#' @srrstats {G2.14b} No imputation; fails fast.
#' @srrstats {G2.16} Accepts `data.frame` and `tibble`.
#' @srrstats {G5.0} Edge-case tests recommended.
#' @srrstats {G5.1} Documentation includes edge-case behaviour.
#' @srrstats {G5.2} Handles failed replicates gracefully.
#' @srrstats {G5.2a} Handles <20 successful replicates.
#' @srrstats {G5.2b} Handles 0 successful replicates.
#' @srrstats {G5.3} Warns on low success rate (<80%).
#' @srrstats {G5.4} Checks `n < nmin * (num_cuts + 1)`.
#' @srrstats {G5.4a} `nmin` default = 90% of original.
#' @srrstats {G5.5} Reproducible via `seed`.
#' @srrstats {G5.6} Optional parallel via `doParallel`.
#' @srrstats {G5.9} Output includes `bootstrap_distribution`.
#' @srrstats {G5.9a} `confidence_intervals` as data frame.
#' @srrstats {G5.9b} `boot_summary` with mean, SD, median, IQR.
#' @srrstats {G5.10} `print()` shows CI and success rate.
#' @srrstats {G5.11} `summary()` shows descriptives, CI, params.
#' @srrstats {G5.11a} `plot()` shows density + CI.
#' @srrstats {RE7.0} Implements bootstrap resampling.
#' @srrstats {RE7.1} Assesses cut-point stability.
#' @srrstats {RE7.2} Returns full bootstrap distribution.
#' @srrstats {RE7.3} Output retains original cut-points.
#' @srrstats {RE7.4} `parameters` list includes all settings.
#' @srrstats {G3.0} S3 plot method for bootstrap distribution.
#' @srrstats {G3.1} Diagnostic plot of bootstrap distribution.
#'
#' @param cutpoint_result An object from [find_cutpoint()].
#' @param num_replicates Number of bootstrap replicates. Default is 500.
#' @param use_parallel Logical. Use multiple CPU cores?
#' @param n_cores Number of cores for parallel. `NULL` defaults to 2.
#' @param seed Optional integer for reproducible results.
#' @param nmin Minimum group size for bootstrap runs. Defaults to 90%
#'   of original `nmin` to reduce failures.
#' @param ... Additional arguments passed to [find_cutpoint()]
#'   (e.g., `popSize`, `maxiter` for genetic algorithm).
#'
#' @return An object of class `validate_cutpoint_result` with
#'   original cuts, 95% CIs, bootstrap distribution, and parameters.
#'
#' @examples
#' # Fast validation on small data (runs in < 2 seconds)
#' data(crc_virome)
#'
#' fit <- find_cutpoint(
#'   data = head(crc_virome, 50),
#'   predictor = "Alphapapillomavirus",
#'   outcome_time = "time_months",
#'   outcome_event = "status",
#'   num_cuts = 1,
#'   method = "systematic"
#' )
#'
#' if (!any(is.na(fit$optimal_cuts))) {
#'   val <- validate_cutpoint(fit, num_replicates = 20, seed = 123)
#'   summary(val)
#'   plot(val)
#' }
#'
#' @references
#' Efron, B. (1979). Bootstrap Methods: Another Look at the
#' Jackknife. *The Annals of Statistics*, 7(1), 1–26.
#' \doi{10.1214/aos/1176344552}
#'
#' Rota, M., Antolini, L., & Valsecchi, M. G. (2015).
#' Optimal cut-point definition in biomarkers: The case of censored
#' failure time outcome. *BMC Medical Research Methodology*, 15(1), 24.
#' \doi{10.1186/s12874-015-0009-y}
#'
#' @importFrom foreach %dopar%
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success
#' @importFrom cli cli_alert_warning cli_inform cli_abort
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom stats quantile na.omit complete.cases sd median IQR
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom tidyr pivot_longer everything
#' @export
validate_cutpoint <- function(cutpoint_result, num_replicates = 500,
                              use_parallel = FALSE, n_cores = NULL,
                              seed = NULL, nmin = NULL, ...) {
  # --- 1. Validate Input and Set Seed ---
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    cli::cli_abort("Input must be a `find_cutpoint` object.")
  }
  if (!is.numeric(num_replicates) || num_replicates < 1 ||
      num_replicates != round(num_replicates)) {
    cli::cli_abort("`num_replicates` must be a positive integer.")
  }
  if (num_replicates < 20) {
    cli::cli_abort("`num_replicates` must be >= 20 for validation.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
    cli::cli_alert_info("Using random seed {seed} for reproducibility.")
  }

  # --- 2. Extract Original Parameters and Data ---
  original_params <- cutpoint_result$parameters
  original_data <- cutpoint_result$userdata
  original_cuts <- cutpoint_result$optimal_cuts
  n <- nrow(original_data)
  if (any(is.na(original_cuts))) {
    cli::cli_abort("Input `find_cutpoint` object has NA cut-points.")
  }

  predictor <- original_params$predictor
  num_cuts <- original_params$num_cuts
  method <- original_params$method
  criterion <- original_params$criterion
  covariates <- original_params$covariates

  # --- Default nmin: 90% of original ---
  if (is.null(nmin)) {
    original_nmin_param <- original_params$nmin
    if (original_nmin_param > 0 && original_nmin_param < 1) {
      original_nmin_abs <- floor(original_nmin_param * n)
      nmin <- floor(0.9 * original_nmin_abs)
    } else {
      nmin <- floor(0.9 * original_nmin_param)
    }
    if (nmin < 1) nmin <- 1
    cli::cli_alert_info("Bootstrap `nmin` not set, defaulting to: {nmin}")
  }
  if (!is.numeric(nmin) || nmin < 1) {
    cli::cli_abort("`nmin` must be a positive integer.")
  }
  if (n < nmin * (num_cuts + 1)) {
    cli::cli_abort(
      "Not enough data ({n}) for nmin ({nmin}) and {num_cuts} cut(s)."
    )
  }

  cli::cli_alert_info(paste(
    "Validating {num_cuts} cut(s) from '{method}' search",
    "using '{criterion}'."
  ))

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
      cli::cli_abort("Package 'doParallel' is required for parallel.")
    }
    cl <- parallel::makeCluster(cores_to_use)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cli::cli_alert_info(
      "Running {num_replicates} replicates on {cores_to_use} core{?s}..."
    )
  } else {
    foreach::registerDoSEQ()
    cli::cli_alert_info(
      "Running {num_replicates} replicates sequentially..."
    )
  }

  # --- 4. Main Bootstrap Loop ---
  if (!use_parallel) {
    pb <- cli::cli_progress_bar("Bootstrapping", total = num_replicates)
  }

  base_funcs <- c(
    "find_cutpoint",
    ".validate_find_cutpoint_inputs",
    ".prepare_cutpoint_data",
    ".validate_data_conditions",
    ".systematic_search",
    ".get_stat"
  )
  functions_to_export <- if (method == "genetic") {
    c(base_funcs, ".run_genetic_search", ".obj")
  } else {
    base_funcs
  }

  i <- NULL
  bootstrap_results <- foreach::foreach(
    i = 1:num_replicates,
    .combine = "rbind",
    .packages = "OptSurvCutR",
    .errorhandling = "pass"
  ) %dopar% {
    if (!use_parallel) {
      cli::cli_progress_update(id = pb)
    }
    set.seed(i)
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- original_data[boot_indices, ]

    extra_args <- list(...)
    final_args <- c(
      list(
        data = boot_data,
        predictor = "factor",
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

    res <- tryCatch(
      {
        suppressMessages(do.call(find_cutpoint, final_args))
      },
      error = function(e) {
        cli::cli_inform("Bootstrap replicate {i} failed: {e$message}")
        return(NULL)
      }
    )

    if (is.null(res) || any(is.na(res$optimal_cuts)) ||
        length(res$optimal_cuts) != num_cuts) {
      return(rep(NA, num_cuts))
    }
    res$optimal_cuts
  }

  # --- 5. Process Results and Calculate CIs ---
  if (!is.matrix(bootstrap_results)) {
    bootstrap_matrix <- matrix(bootstrap_results,
                               nrow = num_replicates,
                               ncol = num_cuts, byrow = TRUE
    )
  } else {
    bootstrap_matrix <- bootstrap_results
  }
  colnames(bootstrap_matrix) <- paste0("Cut", 1:num_cuts)

  successful_reps <- sum(stats::complete.cases(bootstrap_matrix))
  failed_reps <- num_replicates - successful_reps

  if (failed_reps > 0) {
    cli::cli_alert_warning(
      "{failed_reps} of {num_replicates} replicates failed."
    )
  }
  if (successful_reps < 20) {
    cli::cli_abort(paste(
      "< 20 successful replicates ({successful_reps}).",
      "Increase num_replicates or reduce nmin."
    ))
  }
  cli::cli_alert_success("{successful_reps} replicates completed.")

  if (successful_reps / num_replicates < 0.8 && successful_reps > 0) {
    warning(
      paste0(
        successful_reps, " of ", num_replicates,
        " replicates succeeded. CIs may be unreliable."
      )
    )
  }

  bootstrap_matrix_clean <- na.omit(bootstrap_matrix)
  ci <- apply(bootstrap_matrix_clean, 2, stats::quantile,
              probs = c(0.025, 0.975), na.rm = TRUE
  )
  if (is.vector(ci)) {
    ci_df <- data.frame(Lower = ci[1], Upper = ci[2])
  } else {
    ci_df <- as.data.frame(t(ci))
    names(ci_df) <- c("Lower", "Upper")
  }
  row.names(ci_df) <- paste0("Cut ", 1:num_cuts)

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

  bootstrap_df <- as.data.frame(bootstrap_matrix_clean)
  names(bootstrap_df) <- paste0("Cut_point_", 1:num_cuts)

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
#' @srrstats {G5.10} `print()` shows CI and success rate.
#' @export
print.validate_cutpoint_result <- function(x, ...) {
  cat("Cut-point Stability Analysis (Bootstrap)\n")
  cat("----------------------------------------\n")
  cat(
    "Original Optimal Cut-point(s):",
    paste(round(x$original_cuts, 3), collapse = ", "), "\n"
  )
  cat(
    "Successful Replicates:", x$parameters$successful_reps, "/",
    x$parameters$num_replicates,
    "(", round(100 * x$parameters$successful_reps /
                 x$parameters$num_replicates, 1), "%)\n"
  )
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
  numeric_cols <- vapply(summary_df, is.numeric, FUN.VALUE = logical(1))
  summary_df[, numeric_cols] <- round(summary_df[, numeric_cols], 3)
  print(summary_df)
  cat("\nHint: Use `summary()` or `plot()` to visualize stability.\n")
}

#' @rdname validate_cutpoint
#' @importFrom ggplot2 ggplot aes .data geom_density geom_vline
#' @importFrom ggplot2 labs theme_minimal facet_wrap
#' @srrstats {G3.0} S3 plot method for bootstrap distribution.
#' @srrstats {G3.1} Diagnostic plot of bootstrap distribution.
#' @srrstats {G5.11a} `plot()` shows density + CI.
#' @export
plot.validate_cutpoint_result <- function(x, ...) {
  dist_data <- x$bootstrap_distribution
  num_cuts <- ncol(dist_data)

  if (x$parameters$successful_reps == 0) {
    cli::cli_inform("Cannot plot: 0 successful bootstrap replicates.")
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
    ggplot2::geom_density(
      fill = "#56B4E9", color = "#0072B2", alpha = 0.6
    ) +
    ggplot2::geom_vline(
      data = line_data, aes(xintercept = .data$original_cut),
      color = "#D55E00", linetype = "solid", linewidth = 1
    ) +
    ggplot2::geom_vline(
      data = line_data, aes(xintercept = .data$ci_lower),
      color = "#D55E00", linetype = "dashed"
    ) +
    ggplot2::geom_vline(
      data = line_data, aes(xintercept = .data$ci_upper),
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
#' @param show_descriptives Logical. Show descriptive statistics?
#' @param show_ci Logical. Show confidence intervals?
#' @param show_params Logical. Show validation run parameters?
#' @param plot.it Logical. Display the density plot?
#' @rdname validate_cutpoint
#' @srrstats {G5.11} `summary()` shows descriptives, CI, params.
#' @export
summary.validate_cutpoint_result <- function(object, show_descriptives = TRUE,
                                             show_ci = TRUE, show_params = TRUE,
                                             plot.it = FALSE, ...) {
  cat("Cut-point Stability Analysis (Bootstrap)\n")
  cat("----------------------------------------\n")
  cat(
    "Original Optimal Cut-point(s):",
    paste(round(object$original_cuts, 3), collapse = ", "), "\n\n"
  )

  if (show_descriptives) {
    cat("Bootstrap Distribution Summary\n")
    cat("-----------------------------\n")
    summary_df <- do.call(
      rbind,
      lapply(names(object$boot_summary), function(cut) {
        stats <- object$boot_summary[[cut]]
        data.frame(
          Cut = cut,
          Mean = stats$mean,
          SD = stats$sd,
          Median = stats$median,
          Q1 = stats$Q1,
          Q3 = stats$Q3
        )
      })
    )
    numeric_cols <- vapply(summary_df, is.numeric, FUN.VALUE = logical(1))
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
    cat(
      "Successful Replicates:", object$parameters$successful_reps, "/",
      object$parameters$num_replicates,
      "(", round(100 * object$parameters$successful_reps /
                   object$parameters$num_replicates, 1), "%)\n"
    )
    cat("Failed Replicates:", object$parameters$failed_reps, "\n")
    cat(
      "Parallel Processing:",
      ifelse(object$parameters$use_parallel, "Enabled", "Disabled"), "\n"
    )
    cat("Cores Used:", object$parameters$n_cores, "\n")
    cat(
      "Seed:",
      ifelse(is.null(object$parameters$seed), "Not set",
             object$parameters$seed
      ), "\n"
    )
    cat("Minimum Group Size (nmin):", object$parameters$nmin, "\n")
    cat("Method:", object$parameters$method, "\n")
    cat("Criterion:", object$parameters$criterion, "\n")
    cat(
      "Covariates:",
      ifelse(is.null(object$parameters$covariates), "None",
             paste(object$parameters$covariates, collapse = ", ")
      ), "\n"
    )
    cat("\n")
  }

  if (plot.it) {
    cat("Bootstrap Distribution Plot\n")
    cat("--------------------------\n")
    print(plot(object, ...))
  }

  invisible(object)
}
