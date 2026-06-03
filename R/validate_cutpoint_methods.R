# ===================================================================
# S3 METHODS
# Print, Summary, and Plot methods for all 'validate_cutpoint' objects
# ===================================================================

#' @rdname validate_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.17} `print()` shows CI and success rate.
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
    data.frame(Cut = cut, Mean = stats$mean, SD = stats$sd, Median = stats$median, Q1 = stats$Q1, Q3 = stats$Q3)
  }))
  numeric_cols <- vapply(summary_df, is.numeric, FUN.VALUE = logical(1))
  summary_df[, numeric_cols] <- round(summary_df[, numeric_cols], 3)
  print(summary_df)
  cat("\nHint: Use `summary()` or `plot()` to visualise stability.\n")
}

#' @rdname validate_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE6.0} Plot method provided for bootstrap distribution.
#' @export
plot.validate_cutpoint_result <- function(x, ...) {
  dist_data <- x$bootstrap_distribution
  num_cuts <- ncol(dist_data)
  if (x$parameters$successful_reps == 0) {
    cli::cli_inform("Cannot plot: 0 successful bootstrap replicates.")
    return(invisible(NULL))
  }
  plot_data <- tidyr::pivot_longer(dist_data, cols = tidyr::everything(), names_to = "Cut", values_to = "Value", names_prefix = "Cut_point_")
  plot_data$Cut <- paste("Cut-point", plot_data$Cut)
  line_data <- data.frame(
    Cut = paste("Cut-point", 1:num_cuts),
    original_cut = x$original_cuts,
    ci_lower = x$confidence_intervals$Lower,
    ci_upper = x$confidence_intervals$Upper
  )
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Value)) +
    ggplot2::geom_density(fill = "#56B4E9", color = "#0072B2", alpha = 0.6) +
    ggplot2::geom_vline(data = line_data, ggplot2::aes(xintercept = .data$original_cut), color = "#D55E00", linetype = "solid", linewidth = 1) +
    ggplot2::geom_vline(data = line_data, ggplot2::aes(xintercept = .data$ci_lower), color = "#D55E00", linetype = "dashed") +
    ggplot2::geom_vline(data = line_data, ggplot2::aes(xintercept = .data$ci_upper), color = "#D55E00", linetype = "dashed") +
    ggplot2::labs(title = "Bootstrap Distribution of Optimal Cut-points",
                  subtitle = paste(x$parameters$successful_reps, "successful replicates"),
                  x = "Cut-point Value", y = "Density") +
    ggplot2::theme_minimal(base_size = 14)
  if (num_cuts > 1) {
    p <- p + ggplot2::facet_wrap(~Cut, scales = "free_x")
  }
  return(p)
}

#' @rdname validate_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.18} `summary()` shows descriptives, CI, params, and potential over-fitting warnings.
#' @export
summary.validate_cutpoint_result <- function(object, show_descriptives = TRUE, show_ci = TRUE, show_params = TRUE, plot.it = FALSE, ...) {
  cat("Cut-point Stability Analysis (Bootstrap)\n")
  cat("----------------------------------------\n")
  cat("Original Optimal Cut-point(s):", paste(round(object$original_cuts, 3), collapse = ", "), "\n\n")

  if (show_descriptives) {
    cat("Bootstrap Distribution Summary\n-----------------------------\n")
    summary_df <- do.call(rbind, lapply(names(object$boot_summary), function(cut) {
      stats <- object$boot_summary[[cut]]
      data.frame(Cut = cut, Mean = stats$mean, SD = stats$sd, Median = stats$median, Q1 = stats$Q1, Q3 = stats$Q3)
    }))
    numeric_cols <- vapply(summary_df, is.numeric, FUN.VALUE = logical(1))
    summary_df[, numeric_cols] <- round(summary_df[, numeric_cols], 3)
    print(summary_df)
    cat("\n")
  }

  if (show_ci) {
    cat("95% Confidence Intervals\n------------------------\n")
    print(round(object$confidence_intervals, 3))
    cat("\n")
  }

  if (show_params) {
    cat("Validation Parameters\n---------------------\n")
    cat("Replicates Requested:", object$parameters$num_replicates, "\n")
    cat("Successful Replicates:", object$parameters$successful_reps, "/", object$parameters$num_replicates,
        "(", round(100 * object$parameters$successful_reps / object$parameters$num_replicates, 1), "%)\n")
    cat("Failed Replicates:", object$parameters$failed_reps, "\n")
    cat("Cores Used:", object$parameters$n_cores, "\n")
    cat("Seed:", ifelse(is.null(object$parameters$seed), "Not set", object$parameters$seed), "\n")
    cat("Minimum Group Size (nmin):", object$parameters$nmin, "\n")
    cat("Method:", object$parameters$method, "\n")
    cat("Criterion:", object$parameters$criterion, "\n")
    cat("Covariates:", ifelse(is.null(object$parameters$covariates), "None", paste(object$parameters$covariates, collapse = ", ")), "\n\n")
  }

  if (plot.it) {
    cat("Bootstrap Distribution Plot\n--------------------------\n")
    print(plot(object, ...))
  }

  # ===================================================================
  # 3-TIER STABILITY ASSESSMENT (10th-90th Percentile & Override)
  # ===================================================================

  # 1. Safely calculate the 10th-90th Percentile spread
  if (!is.null(object$userdata) && !is.null(object$userdata$factor)) {
    p10 <- stats::quantile(object$userdata$factor, 0.10, na.rm = TRUE)
    p90 <- stats::quantile(object$userdata$factor, 0.90, na.rm = TRUE)
    data_spread <- p90 - p10
  } else if (!is.null(object$bootstrap_distribution)) {
    p10 <- stats::quantile(object$bootstrap_distribution, 0.10, na.rm = TRUE)
    p90 <- stats::quantile(object$bootstrap_distribution, 0.90, na.rm = TRUE)
    data_spread <- p90 - p10
  } else {
    data_spread <- NA
  }

  # 2. Extract medians
  cut_names <- rownames(object$confidence_intervals)
  medians <- vapply(cut_names, function(n) {
    clean_n <- gsub("[^0-9]", "", n)
    match_idx <- grep(clean_n, names(object$boot_summary))
    if(length(match_idx) == 0) return(NA)
    object$boot_summary[[match_idx]]$median
  }, FUN.VALUE = numeric(1))

  # Failsafe: If the spread is 0, missing, or Infinite, fall back to the magnitude of the medians
  if (is.na(data_spread) || is.infinite(data_spread) || data_spread <= 0) {
    data_spread <- suppressWarnings(max(abs(medians), na.rm = TRUE))
    if (data_spread == 0 || is.infinite(data_spread)) data_spread <- 1
  }

  # 3. Calculate relative widths using the safe 80% spread
  lower_ci <- object$confidence_intervals$Lower
  upper_ci <- object$confidence_intervals$Upper

  relative_widths <- (upper_ci - lower_ci) / data_spread
  max_rciw <- suppressWarnings(max(relative_widths, na.rm = TRUE))

  # Prevent multiplying by NA/Inf
  if (is.finite(max_rciw)) {
    stability_pct <- round(max_rciw * 100, 1)
  } else {
    stability_pct <- NA
  }

  # Safely extract the worst cut name
  worst_cut_name <- "Unknown Cut"
  valid_idx <- which.max(relative_widths)
  if (length(valid_idx) > 0) {
    worst_cut_name <- rownames(object$confidence_intervals)[valid_idx]
  }

  # 4. Check for Biological Override (Non-overlapping CIs)
  is_perfectly_separated <- FALSE
  if (nrow(object$confidence_intervals) > 1) {
    is_perfectly_separated <- TRUE

    # Safe sort: Order by the lower bound to ensure correct [i] vs [i+1] comparison
    ci_ordered <- object$confidence_intervals[order(object$confidence_intervals$Lower), ]

    for (i in seq_len(nrow(ci_ordered) - 1)) {
      # Safely handle NAs in the overlap check
      if (is.na(ci_ordered$Upper[i]) || is.na(ci_ordered$Lower[i + 1])) {
        is_perfectly_separated <- FALSE
        break
      }
      if (ci_ordered$Upper[i] >= ci_ordered$Lower[i + 1]) {
        is_perfectly_separated <- FALSE
        break
      }
    }
  }

  # 5. Print the Final Assessment
  cat("\nStability Assessment:\n")
  cat("---------------------\n")

  if (!is.finite(max_rciw) || is.na(stability_pct)) {
    cli::cli_alert_warning("Stability could not be calculated (missing or infinite CI data).")

  } else {
    cli::cli_text("Maximum CI Width (Relative to 10th-90th Percentile Range): {.strong {stability_pct}%}\n")

    # Determine structural properties for the 4-Tier Logic
    is_multi_cut <- nrow(object$confidence_intervals) > 1
    has_overlap <- is_multi_cut && !is_perfectly_separated
    has_distinct_separation <- is_multi_cut && is_perfectly_separated

    if (max_rciw < 0.30) {
      if (has_overlap) {
        # TIER 3: OVERLAP DOWNGRADE (Extremely tight variance, but they still overlap)
        cli::cli_alert_warning("Model Status: CAUTION (Tier 3) - OVERLAP DOWNGRADE")
        cli::cli_text("The mathematical variance is very low ({stability_pct}%), but the Confidence Intervals overlap.")
        cli::cli_text("This indicates that while the thresholds are mathematically stable, the resulting risk cohorts blend together in the 'grey zones'.")
        cli::cli_bullets(c(
          "*" = "If distinct separation is required for decision-making, consider reducing {.arg num_cuts}.",
          "*" = "If exploratory, this model is acceptable but should be interpreted with caution in the overlapping ranges."
        ))
      } else {
        # TRUE TIER 1: OPTIMAL
        cli::cli_alert_success("Model Status: OPTIMAL (Tier 1)")
        cli::cli_text("The thresholds are highly consistent across samples with clean separation between risk cohorts.")
      }

    } else if (max_rciw >= 0.30 && max_rciw <= 0.60) {
      if (has_distinct_separation) {
        # TIER 2: DISTINCT (Moderate Variance, Zero Overlap)
        cli::cli_alert_success("Model Status: DISTINCT (Tier 2)")
        cli::cli_text("The relative mathematical variance is moderate ({stability_pct}%), but the 95% Confidence Intervals for your cut-points do not overlap.")
        cli::cli_text("This indicates the algorithm found mathematically distinct subpopulations despite exact threshold variance.")
      } else if (has_overlap) {
        # TIER 3: CAUTION (Moderate Variance, Overlapping)
        cli::cli_alert_warning("Model Status: CAUTION (Tier 3)")
        cli::cli_text("Moderate instability detected ({stability_pct}%), and Confidence Intervals overlap.")
        cli::cli_text("This indicates that the resulting risk cohorts blend together in the 'grey zones'.")
        cli::cli_bullets(c(
          "*" = "If distinct separation is required for decision-making, consider reducing {.arg num_cuts}.",
          "*" = "If exploratory, this model is acceptable but should be interpreted with caution in the overlapping ranges."
        ))
      } else {
        # TIER 3: CAUTION (Single Cut Model)
        cli::cli_alert_warning("Model Status: CAUTION (Tier 3)")
        cli::cli_text("Moderate instability detected ({stability_pct}%). The threshold is sensitive to sample variance but remains within an acceptable range.")
      }

    } else {
      # max_rciw > 0.60
      if (has_distinct_separation) {
        # TIER 2: DISTINCT (High Variance Override, Zero Overlap)
        cli::cli_alert_success("Model Status: DISTINCT (Tier 2) - SEPARATION OVERRIDE")
        cli::cli_text("The relative mathematical variance is high ({stability_pct}%), but the 95% Confidence Intervals for your cut-points do not overlap.")
        cli::cli_text("This indicates the algorithm found mathematically distinct, highly stable subpopulations despite a narrow data range.")
      } else {
        # TIER 4: UNSTABLE (High Variance, Overlapping or Single Cut)
        cli::cli_alert_danger("Model Status: UNSTABLE (Tier 4)")
        if (has_overlap) {
          cli::cli_text("High instability detected ({stability_pct}%)! The cut-points are highly sensitive to sample changes and overlap significantly, indicating potential over-fitting to noise.")
        } else {
          cli::cli_text("High instability detected ({stability_pct}%)! The cut-point is highly sensitive to sample changes, indicating potential over-fitting to noise.")
        }
        cli::cli_bullets(c(
          "!" = "The primary source of instability is {.strong {worst_cut_name}}.",
          "x" = "Recommendation: Reduce {.arg num_cuts} or increase {.arg nmin}.",
          "->" = "See the Rescue Protocol: {.code vignette('troubleshooting', package = 'OptSurvCutR')}"
        ))
      }
    }
  }
  cat("\n")
}
