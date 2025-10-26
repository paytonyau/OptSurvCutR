#' Plot Optimization Curve from a Systematic Search
#'
#' @description
#' Visualizes the optimization process from a `find_cutpoint()` systematic search
#' with `num_cuts = 1`. It plots the chosen metric (e.g., Log-Rank statistic,
#' Hazard Ratio, or p-value) against all evaluated cut-points.
#'
#' This plot helps to visually confirm the optimal cut-point and assess the
#' sensitivity of the metric to the choice of cut-point.
#'
#' @param cutpoint_result An object returned by `find_cutpoint()` where
#'   `method = "systematic"` and `num_cuts = 1` was used.
#'
#' @return A ggplot object showing the metric values across the range of
#'   cut-points. The optimal cut-point is marked with a vertical dashed line.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline labs theme_minimal .data
#' @importFrom cli cli_abort
#' @importFrom tools toTitleCase
#' @export
plot_optimization_curve <- function(cutpoint_result) {

  # --- 1. Input Validation ---
  params <- cutpoint_result$parameters
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    cli::cli_abort("Input must be an object from the {.fn find_cutpoint} function.")
  }
  if (is.null(params$method) || params$method != "systematic") {
    cli::cli_abort("This plot is only for results from {.fn find_cutpoint} with `method = \"systematic\"`.")
  }
  if (is.null(params$num_cuts) || params$num_cuts != 1) {
    cli::cli_abort("This plot is only supported for results with `num_cuts = 1`.")
  }
  if (is.null(cutpoint_result$all_stats) || nrow(cutpoint_result$all_stats) == 0) {
    cli::cli_abort("The {.arg cutpoint_result} object does not contain the necessary `all_stats` data for this plot.")
  }

  # --- 2. Data Preparation ---
  plot_data <- cutpoint_result$all_stats
  optimal_cut <- cutpoint_result$optimal_cuts[1]
  criterion <- params$criterion

  # Determine plot labels based on the criterion used
  y_label <- switch(criterion,
                    "logrank" = "Log-Rank Statistic",
                    "hazard_ratio" = "Hazard Ratio (HR)",
                    "p_value" = "P-value",
                    tools::toTitleCase(criterion)
  )
  plot_title <- paste(y_label, "vs. Cut-point")
  subtitle_text <- if (!is.na(optimal_cut)) {
    paste("Optimal cut-point at", round(optimal_cut, 3))
  } else {
    "No optimal cut-point found."
  }

  # --- 3. Generate Plot ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$cut1, y = .data$stat)) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
    ggplot2::labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "Cut-point Value",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 14)

  # Add line for optimal cut if it was found
  if (!is.na(optimal_cut)) {
    p <- p + ggplot2::geom_vline(
      xintercept = optimal_cut,
      linetype = "dashed",
      color = "#D55E00",
      linewidth = 1.2
    )
  }

  # For HR plots, add a reference line at 1
  if (criterion == "hazard_ratio") {
    p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dotted")
  }

  return(p)
}
