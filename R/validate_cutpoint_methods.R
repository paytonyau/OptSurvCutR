# ===================================================================
# S3 METHODS
# Print, Summary, and Plot methods for all 'validate_cutpoint' objects
# ===================================================================

#' @rdname plot_validation
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.17} \code{print()} shows CI and success rate.
#' @export
print.validate_cutpoint_result <- function(x, ...) {
  cat("Cut-point Stability Analysis (Bootstrap)\n")
  cat("----------------------------------------\n")
  cat("Original Optimal Cut-point(s):", paste(round(x$original_cuts, 3), collapse = ", "), "\n")
  cat(
    "Successful Replicates:", x$parameters$successful_reps, "/", x$parameters$num_replicates,
    "(", round(100 * x$parameters$successful_reps / x$parameters$num_replicates, 1), "%)\n"
  )
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
  invisible(x)
}

#' @rdname plot_validation
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
  
  plot_data <- tidyr::pivot_longer(dist_data, cols = tidyr::everything(), 
                                   names_to = "Cut", values_to = "Value", 
                                   names_prefix = "Cut_point_")
  plot_data$Cut <- paste("Cut-point", plot_data$Cut)
  
  line_data <- data.frame(
    Cut = paste("Cut-point", seq_len(num_cuts)),
    original_cut = x$original_cuts,
    ci_lower = x$confidence_intervals$Lower,
    ci_upper = x$confidence_intervals$Upper
  )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$Value)) +
    ggplot2::geom_density(fill = "#56B4E9", color = "#0072B2", alpha = 0.6) +
    ggplot2::geom_vline(data = line_data, ggplot2::aes(xintercept = .data$original_cut), 
                        color = "#D55E00", linetype = "solid", linewidth = 1) +
    ggplot2::geom_vline(data = line_data, ggplot2::aes(xintercept = .data$ci_lower), 
                        color = "#D55E00", linetype = "dashed") +
    ggplot2::geom_vline(data = line_data, ggplot2::aes(xintercept = .data$ci_upper), 
                        color = "#D55E00", linetype = "dashed") +
    ggplot2::labs(
      title = "Bootstrap Distribution of Optimal Cut-points",
      subtitle = paste(x$parameters$successful_reps, "successful replicates"),
      x = "Cut-point Value", y = "Density"
    ) +
    ggplot2::theme_minimal(base_size = 14)
  
  if (num_cuts > 1) {
    p <- p + ggplot2::facet_wrap(~Cut, scales = "free_x")
  }
  return(p)
}

#' @rdname plot_validation
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.18} `summary()` shows descriptives, CI, params, and potential over-fitting warnings.
#' @export
summary.validate_cutpoint_result <- function(object, show_descriptives = TRUE, show_ci = TRUE, 
                                             show_params = TRUE, plot.it = FALSE, ...) {
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
    cat(
      "Successful Replicates:", object$parameters$successful_reps, "/", object$parameters$num_replicates,
      "(", round(100 * object$parameters$successful_reps / object$parameters$num_replicates, 1), "%)\n"
    )
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
  # 4-TIER STABILITY ASSESSMENT ENGINE
  # ===================================================================
  
  if (!is.null(object$userdata) && !is.null(object$userdata$factor)) {
    p10 <- stats::quantile(object$userdata$factor, 0.10, na.rm = TRUE)
    p90 <- stats::quantile(object$userdata$factor, 0.90, na.rm = TRUE)
    data_spread <- p90 - p10
  } else if (!is.null(object$bootstrap_distribution)) {
    # ✅ FIXED: Force drop vector evaluation to secure quantile parsing across tibble structures
    target_vec <- object$bootstrap_distribution[, 1, drop = TRUE]
    p10 <- stats::quantile(target_vec, 0.10, na.rm = TRUE)
    p90 <- stats::quantile(target_vec, 0.90, na.rm = TRUE)
    data_spread <- p90 - p10
  } else {
    data_spread <- NA_real_
  }
  
  cut_names <- rownames(object$confidence_intervals)
  medians <- vapply(cut_names, function(n) {
    clean_n <- gsub("[^0-9]", "", n)
    # ✅ FIXED: Exact name reconstruction prevents grep collisions on high-dimensional models (e.g. Cut 1 vs Cut 10)
    target_name <- paste0("Cut", clean_n)
    match_idx <- which(names(object$boot_summary) == target_name)
    if (length(match_idx) == 0) {
      return(NA_real_)
    }
    object$boot_summary[[match_idx]]$median
  }, FUN.VALUE = numeric(1))
  
  if (is.na(data_spread) || is.infinite(data_spread) || data_spread <= 0) {
    data_spread <- suppressWarnings(max(abs(medians), na.rm = TRUE))
    if (data_spread == 0 || is.infinite(data_spread)) data_spread <- 1
  }
  
  lower_ci <- object$confidence_intervals$Lower
  upper_ci <- object$confidence_intervals$Upper
  
  relative_widths <- (upper_ci - lower_ci) / data_spread
  max_rciw <- suppressWarnings(max(relative_widths, na.rm = TRUE))
  
  if (is.finite(max_rciw)) {
    stability_pct <- round(max_rciw * 100, 1)
  } else {
    stability_pct <- NA_real_
  }
  
  worst_cut_name <- "Unknown Cut"
  valid_idx <- which.max(relative_widths)
  if (length(valid_idx) > 0) {
    worst_cut_name <- rownames(object$confidence_intervals)[valid_idx]
  }
  
  is_perfectly_separated <- FALSE
  if (nrow(object$confidence_intervals) > 1) {
    is_perfectly_separated <- TRUE
    ci_ordered <- object$confidence_intervals[order(object$confidence_intervals$Lower), ]
    
    for (i in seq_len(nrow(ci_ordered) - 1)) {
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
  
  cat("\nStability Assessment:\n")
  cat("---------------------\n")
  
  if (!is.finite(max_rciw) || is.na(stability_pct)) {
    cli::cli_alert_warning("Stability could not be calculated (missing or infinite CI data).")
  } else {
    cli::cli_text("Maximum CI Width (Relative to 10th-90th Percentile Range): {.strong {stability_pct}%}\n")
    
    is_multi_cut <- nrow(object$confidence_intervals) > 1
    has_overlap <- is_multi_cut && !is_perfectly_separated
    has_distinct_separation <- is_multi_cut && is_perfectly_separated
    
    if (max_rciw < 0.30) {
      if (has_overlap) {
        cli::cli_alert_warning("Model Status: CAUTION (Tier 3) - OVERLAP DOWNGRADE")
        cli::cli_text("The mathematical variance is very low ({stability_pct}%), but the Confidence Intervals overlap.")
        cli::cli_bullets(c(
          "*" = "Consider reducing {.arg num_cuts} if distinct separation is required for clinical application."
        ))
      } else {
        cli::cli_alert_success("Model Status: OPTIMAL (Tier 1)")
        cli::cli_text("The thresholds are highly consistent across samples with clean separation between risk cohorts.")
      }
    } else if (max_rciw >= 0.30 && max_rciw <= 0.60) {
      if (has_distinct_separation) {
        cli::cli_alert_success("Model Status: DISTINCT (Tier 2)")
        cli::cli_text("The relative mathematical variance is moderate ({stability_pct}%), but 95% Confidence Intervals do not overlap.")
      } else {
        cli::cli_alert_warning("Model Status: CAUTION (Tier 3)")
        cli::cli_text("Moderate instability detected ({stability_pct}%), and Confidence Intervals overlap.")
      }
    } else {
      if (has_distinct_separation) {
        cli::cli_alert_success("Model Status: DISTINCT (Tier 2) - SEPARATION OVERRIDE")
        cli::cli_text("The relative mathematical variance is high ({stability_pct}%), but 95% Confidence Intervals do not overlap.")
      } else {
        cli::cli_alert_danger("Model Status: UNSTABLE (Tier 4)")
        cli::cli_bullets(c(
          "!" = "The primary source of instability is {.strong {worst_cut_name}}.",
          "x" = "Recommendation: Reduce {.arg num_cuts} or increase {.arg nmin}."
        ))
      }
    }
  }
  cat("\n")
  invisible(object)
}