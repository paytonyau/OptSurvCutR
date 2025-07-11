#' Plot Effect Size vs. Cut-point
#'
#' Visualizes how the effect size (Hazard Ratio or Odds Ratio) and its
#' confidence interval change across the full range of possible cut-points.
#'
#' @param cutpoint_result An object returned by `find_cutpoint(method = "systematic")`.
#' @return A ggplot object.
#' @export
plot_effect_size <- function(cutpoint_result) {

  if (!inherits(cutpoint_result, "find_cutpoint_systematic")) {
    stop("This plot is only for results from find_cutpoint(method = 'systematic').", call. = FALSE)
  }

  plot_data <- cutpoint_result$allcut
  optimal_cut <- cutpoint_result$best_by_vote
  analysis_type <- cutpoint_result$parameters$analysis_type

  if (analysis_type == "survival") {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Cut1, y = HR)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = HR_low, ymax = HR_up), alpha = 0.2, fill = "dodgerblue") +
      ggplot2::geom_line(color = "dodgerblue", linewidth = 1) +
      ggplot2::labs(y = "Hazard Ratio (HR)", title = "Hazard Ratio vs. Cut-point")
  } else { # logistic
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Cut1, y = OR)) +
      ggplot2::geom_ribbon(ggplot2::aes(ymin = OR_low, ymax = OR_up), alpha = 0.2, fill = "darkorange") +
      ggplot2::geom_line(color = "darkorange", linewidth = 1) +
      ggplot2::labs(y = "Odds Ratio (OR)", title = "Odds Ratio vs. Cut-point")
  }

  p <- p +
    ggplot2::geom_hline(yintercept = 1, linetype = "dotted") +
    ggplot2::geom_vline(xintercept = optimal_cut, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::labs(
      subtitle = paste("Optimal cut-point at", round(optimal_cut, 2)),
      x = "Cut-point Value"
    ) +
    ggplot2::theme_minimal()

  return(p)
}


#' Plot a Waterfall Chart of Classification
#'
#' Visualizes how an optimal cut-point classifies individual subjects in a
#' binary outcome analysis.
#'
#' @param cutpoint_result An object returned by `find_cutpoint(method = "systematic")`
#'   from a binary outcome analysis.
#' @return A ggplot object.
#' @export
plot_waterfall <- function(cutpoint_result) {

  if (!inherits(cutpoint_result, "find_cutpoint_systematic") ||
      cutpoint_result$parameters$analysis_type != "logistic") {
    stop("Waterfall plot is only for binary outcome results from method = 'systematic'.", call. = FALSE)
  }

  plot_data <- cutpoint_result$userdata
  optimal_cut <- cutpoint_result$best_by_vote

  # Prepare data for plotting
  plot_data <- plot_data %>%
    dplyr::arrange(factor) %>%
    dplyr::mutate(
      patient_id = dplyr::row_number(),
      # Classify based on the cut-point
      classified_group = ifelse(factor <= optimal_cut, 0, 1),
      # Check if classification was correct
      is_correct = (classified_group == outcome)
    )

  # Plot the waterfall chart
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = patient_id, y = factor, fill = is_correct)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_hline(yintercept = optimal_cut, linetype = "dashed", color = "red", linewidth = 1) +
    ggplot2::scale_fill_manual(
      name = "Classification",
      values = c("TRUE" = "forestgreen", "FALSE" = "firebrick"),
      labels = c("TRUE" = "Correct", "FALSE" = "Incorrect")
    ) +
    ggplot2::labs(
      title = "Waterfall Plot of Patient Classification",
      subtitle = paste("Cut-point at", round(optimal_cut, 2)),
      x = "Patients (ordered by predictor value)",
      y = "Predictor Value"
    ) +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())

  return(p)
}
