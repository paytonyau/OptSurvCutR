# ===================================================================
# PLOTTING FUNCTIONS
# Diagnostic and publication-ready plots.
# ===================================================================
#' Plot Optimisation Curve from a Systematic Search
#'
#' @description
#' Plots the metric from a `find_cutpoint()` systematic search
#' (`num_cuts = 1`). It plots the statistic (e.g., Log-Rank, HR,
#' p-value) against all evaluated cut-points.
#'
#' Helps confirm the optimum and assess sensitivity.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.0} Input object class validated with `inherits()`.
#' @srrstats {G2.1} Input types validated.
#' @srrstats {G2.13} Missing/empty `all_stats` triggers `cli_abort()`.
#' @srrstats {G2.14a} `NA` in `optimal_cuts` handled in subtitle.
#' @srrstats {G3.1} Diagnostic plot for systematic search.
#' @srrstats {G5.2} Graceful handling of failed/NA results.
#' @srrstats {G5.2b} Non-fatal messages if plot cannot be produced.
#'
#' @param cutpoint_result A `find_cutpoint` object
#'   (`method = "systematic"`, `num_cuts = 1`).
#'
#' @return A `ggplot` object. The optimal cut-point is marked with a
#'   vertical dashed line.
#'
#' @examples
#' data(crc_virome)
#' fit <- find_cutpoint(
#'   data = head(crc_virome, 40),
#'   predictor = "Alphapapillomavirus",
#'   outcome_time = "time_months",
#'   outcome_event = "status",
#'   num_cuts = 1,
#'   method = "systematic"
#' )
#'
#' if (!any(is.na(fit$optimal_cuts))) {
#'   plot_optimization_curve(fit)
#' }
#'
#' @references
#' Altman, D. G., Lausen, B., Sauerbrei, W., & Schumacher,
#' M. (1994). Dangers of Using "Optimal" Cutpoints in the Evaluation of
#' Prognostic Factors. *JNCI: Journal of the National Cancer Institute*,
#' 86(11), 829-835. \doi{10.1093/jnci/86.11.829}
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline geom_hline
#'   labs theme_minimal .data
#' @importFrom cli cli_abort
#' @importFrom tools toTitleCase
#' @export
plot_optimization_curve <- function(cutpoint_result) {
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    cli::cli_abort(
      "Input must be an object from the {.fn find_cutpoint} function."
    )
  }

  params <- cutpoint_result$parameters
  if (is.null(params)) {
    cli::cli_abort("Missing `parameters` in cutpoint_result.")
  }
  if (params$method != "systematic") {
    cli::cli_abort(
      "This plot is only for `method = \"systematic\"`."
    )
  }
  if (params$num_cuts != 1) {
    cli::cli_abort("This plot is only for `num_cuts = 1`.")
  }

  if (is.null(cutpoint_result$all_stats) ||
    !is.data.frame(cutpoint_result$all_stats)) {
    cli::cli_abort(
      "`cutpoint_result` must have a valid `all_stats` data frame."
    )
  }
  if (nrow(cutpoint_result$all_stats) == 0) {
    cli::cli_abort("`all_stats` is empty; no cuts evaluated.")
  }

  plot_data <- cutpoint_result$all_stats
  optimal_cut <- cutpoint_result$optimal_cuts[1]
  criterion <- params$criterion
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

  p <- ggplot2::ggplot(
    plot_data, ggplot2::aes(x = .data$cut1, y = .data$stat)
  ) +
    ggplot2::geom_line(color = "#0072B2", linewidth = 1) +
    ggplot2::labs(
      title = plot_title,
      subtitle = subtitle_text,
      x = "Cut-point Value",
      y = y_label
    ) +
    ggplot2::theme_minimal(base_size = 14)

  if (!is.na(optimal_cut)) {
    p <- p + ggplot2::geom_vline(
      xintercept = optimal_cut,
      linetype = "dashed",
      color = "#D55E00",
      linewidth = 1.2
    )
  }
  if (criterion == "hazard_ratio") {
    p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dotted")
  }

  return(p)
}

# ===================================================================
#' Cox Proportional Hazards Diagnostic Plot
# ===================================================================
#'
#' @description
#' Plots Schoenfeld residuals to check the proportional hazards (PH)
#' assumption for the final Cox model. Satisfies SRR **G3.1a**.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G3.1a} Provides diagnostic residual plot (Schoenfeld)
#' for Cox PH assumption.
#' @srrstats {G5.2} Graceful handling when no valid cut-points exist.
#' @srrstats {G5.2b} Non-fatal informative message when Cox model fails.
#' @srrstats {G2.0} Input object class validated with `inherits()`.
#' @srrstats {G2.13} Uses `cli_inform()` for non-error user feedback.
#' @srrstats {RE6.1} Includes global Schoenfeld test p-value in plot.
#'
#' @param x A `find_cutpoint` result object.
#' @param ... Additional arguments passed to `survminer::ggcoxzph()`.
#'
#' @return A list of `ggplot` objects (`ggcoxzph` plot).
#'
#' @examples
#' data(crc_virome)
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
#'   plot_cox_diagnostics(fit)
#' }
#'
#' @references
#' Cox, D. R. (1972). Regression Models and Life-Tables. *Journal
#' of the Royal Statistical Society: Series B (Methodological)*, 34(2),
#' 187-202. \doi{10.1111/j.2517-6161.1972.tb00899.x}
#'
#' @importFrom survival coxph cox.zph Surv
#' @importFrom survminer ggcoxzph
#' @importFrom cli cli_inform cli_abort
#' @importFrom stats as.formula
#' @export
plot_cox_diagnostics <- function(x, ...) {
  if (!inherits(x, "find_cutpoint")) {
    cli::cli_abort("Input must be a {.cls find_cutpoint} object.")
  }
  if (is.null(x$optimal_cuts) || any(is.na(x$optimal_cuts))) {
    cli::cli_inform("No valid cut-points; skipping diagnostics.")
    return(invisible(NULL))
  }
  data <- x$userdata
  cuts <- x$optimal_cuts
  data$group <- cut(data$factor,
    breaks = c(-Inf, cuts, Inf),
    labels = paste0("G", 1:(length(cuts) + 1))
  )
  formula_str <- "survival::Surv(time, event) ~ group"
  if (!is.null(x$parameters$covariates)) {
    formula_str <- paste(
      formula_str, "+",
      paste(x$parameters$covariates, collapse = " + ")
    )
  }
  fit <- tryCatch(
    survival::coxph(stats::as.formula(formula_str), data = data),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    cli::cli_inform("Cox model failed; cannot generate diagnostics.")
    return(invisible(NULL))
  }
  zph <- survival::cox.zph(fit)
  p <- survminer::ggcoxzph(zph, ...)
  p <- lapply(p, function(g) {
    g + ggplot2::labs(title = if (is.null(g$labels$title)) {
      "Schoenfeld Residuals"
    } else {
      g$labels$title
    })
  })
  global_p <- round(zph$table[nrow(zph$table), "p"], 4)
  p[[length(p)]] <- p[[length(p)]] +
    ggplot2::labs(subtitle = paste0("Global chi-squared p-value = ", global_p))
  return(p)
}
