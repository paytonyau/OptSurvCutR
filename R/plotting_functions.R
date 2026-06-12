# ===================================================================
# PLOTTING FUNCTIONS
# Diagnostic, publication-ready, and optimized static plots.
# ===================================================================

#' Custom Clinical Theme for OptSurvCutR
#'
#' @description
#' A clean, publication-ready \code{ggplot2} theme used as the standard
#' aesthetic for all OptSurvCutR plots. Features subtle light-grey
#' background gridlines for high-precision tracing.
#'
#' @param base_size Base font size, default is 14.
#' @return A \code{ggplot2} theme object containing customized layout parameters.
#'
#' @importFrom ggplot2 theme_minimal theme element_text element_line rel element_rect
#' @export
#' @examples
#' library(ggplot2)
#' mock_df <- data.frame(x = 1:10, y = 1:10)
#' ggplot(mock_df, aes(x, y)) +
#'     geom_point() +
#'     theme_optsurv()
theme_optsurv <- function(base_size = 14) {
  ggplot2::theme_minimal(base_size = base_size) +
    ggplot2::theme(
      # --- UNIFORM LIGHT-GREY GRIDLINES ---
      panel.grid.major = ggplot2::element_line(color = "grey92", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_line(color = "grey96", linewidth = 0.25),
      
      # --- TYPOGRAPHY & LAYOUT ESTHETICS ---
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.15)),
      plot.subtitle = ggplot2::element_text(color = "grey40", size = ggplot2::rel(0.9)),
      legend.position = "bottom",
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.0)),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(0.95)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.85))
    )
}

#' Master S3 Plot Router for find_cutpoint
#'
#' @description
#' Unified plotting dispatch network routing to clinical survival curves,
#' predictor density distributions, hazard ratio forest charts, multi-dimensional
#' objective surfaces, or conditional landmark stratification assets.
#'
#' @param x A \code{find_cutpoint} result object.
#' @param type Plot framework type: \code{"outcome"}, \code{"distribution"}, \code{"forest"},
#'        \code{"surface"}, \code{"trajectory"}, \code{"diagnostic"}, or \code{"landmark"}.
#' @param return_data Logical. If \code{TRUE}, exits the router early and returns the
#'        assigned underlying data frame template.
#' @param landmark Numeric. The operational milestone timestamp used if \code{type = "landmark"}.
#' @param ... Additional arguments passed down to downstream rendering pipelines.
#' @return A \code{ggplot} canvas object, a multi-panel \code{patchwork} collection, or
#'        a \code{data.frame} if \code{return_data = TRUE}.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE6.0} Implements a robust default plot method for the core model object.
#' @srrstats {RE6.1} Establishes generic S3 plot method routing dispatched directly on model classes.
#'
#' @importFrom cli cli_abort
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   # Build clean local simulation objects with an explicit survival risk split
#'   set.seed(123)
#'   mock_df <- data.frame(
#'      time   = c(runif(15, 50, 100), runif(15, 5, 25)),
#'      event  = rep(1, 30),
#'      factor = c(rnorm(15, 5, 0.5), rnorm(15, 15, 0.5))
#'   )
#'   res <- find_cutpoint(
#'      mock_df, "factor", "time", "event",
#'      num_cuts = 1, method = "systematic", quiet = TRUE, nmin = 3
#'   )
#'   p <- plot(res, type = "distribution")
#' }
plot.find_cutpoint <- function(x, type = c(
  "outcome", "distribution", "forest",
  "surface", "trajectory",
  "diagnostic", "landmark"
),
return_data = FALSE, landmark = NULL, ...) {
  type <- match.arg(type)
  
  if (is.null(x$optimal_cuts) || any(is.na(x$optimal_cuts))) {
    cli::cli_abort("Cannot generate plots: No valid optimal cut-points found in this object.")
  }
  
  df <- x$userdata
  cuts <- sort(x$optimal_cuts)
  num_cuts <- length(cuts)
  
  if (exists("cpp_get_group_assignments", mode = "function") && isTRUE(x$parameters$use_cpp)) {
    grp_vector <- cpp_get_group_assignments(df$factor, cuts)
    df$group <- structure(grp_vector, levels = as.character(1:(num_cuts + 1)), class = "factor")
  } else {
    df$group <- as.factor(findInterval(df$factor, cuts, left.open = TRUE) + 1L)
  }
  
  if (return_data) {
    return(df)
  }
  
  switch(type,
         "outcome"      = .plot_km_curve(x, df, ...),
         "distribution" = .plot_density_cuts(x, df, ...),
         "forest"       = .plot_hr_forest(x, df, ...),
         "surface"      = plot_optimisation_curve(x, ...),
         "trajectory"   = .plot_genetic_trajectory(x, ...),
         "diagnostic"   = plot_cutpoint_residuals(x, ...),
         "landmark"     = plot_landmark_stratification(x, landmark = landmark, ...)
  )
}

#' Plot Optimisation Curve or Surface from Search
#'
#' @description
#' Plots the metric landscape evaluated across coordinates. Maps a 1D optimization line
#' for 1-cut systematic setups, or a 2D topographic grid profile for 2-cut layouts.
#'
#' @param cutpoint_result A \code{find_cutpoint} object.
#' @param ... Unused dots.
#' @return A valid \code{ggplot} object detailing evaluation statistics vs coordinates.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE6.2} Visualizes the continuous fitted values and optimization landscape of the model.
#'
#' @importFrom ggplot2 ggplot aes geom_line geom_vline geom_hline labs geom_tile scale_fill_viridis_c scale_color_manual element_blank
#' @importFrom rlang .data
#' @importFrom cli cli_abort
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   # Build clean simulation objects using a 2-cut systematic search
#'   # to populate the 2D grid matrix log array
#'   set.seed(123)
#'   mock_df <- data.frame(
#'      time   = c(runif(10, 50, 100), runif(10, 20, 60), runif(10, 5, 25)),
#'      event  = rep(1, 30),
#'      factor = c(rnorm(10, 2, 0.2), rnorm(10, 7, 0.2), rnorm(10, 15, 0.2))
#'   )
#'   res <- find_cutpoint(
#'      mock_df, "factor", "time", "event",
#'      num_cuts = 2, method = "systematic", quiet = TRUE, nmin = 3
#'   )
#'   p <- plot_optimisation_curve(res)
#' }
plot_optimisation_curve <- function(cutpoint_result, ...) {
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    cli::cli_abort("Input must be an object from the {.fn find_cutpoint} function.")
  }
  
  params <- cutpoint_result$parameters
  if (params$method != "systematic") {
    cli::cli_abort("Surface mapping (`type = 'surface'`) is only available for `method = 'systematic'`.")
  }
  if (is.null(cutpoint_result$all_stats) || !is.data.frame(cutpoint_result$all_stats)) {
    cli::cli_abort("The results object must contain a valid grid log array in `all_stats`.")
  }
  
  plot_data <- cutpoint_result$all_stats
  criterion <- params$criterion
  
  y_label <- switch(criterion,
                    "logrank" = "Log-Rank Statistic",
                    "hazard_ratio" = "Hazard Ratio",
                    "p_value" = "P-value",
                    criterion
  )
  
  if (params$num_cuts == 1) {
    optimal_cut <- cutpoint_result$optimal_cuts[1]
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$cut1, y = .data$stat)) +
      ggplot2::geom_line(color = "#0072B2", linewidth = 1.2) +
      ggplot2::labs(
        title = paste(y_label, "vs. Cut-point Coordinate"),
        subtitle = paste("Optimal cut-point discovered at:", round(optimal_cut, 3)),
        x = "Cut-point Threshold Location", y = y_label
      ) +
      theme_optsurv()
    
    if (!is.na(optimal_cut)) {
      p <- p + ggplot2::geom_vline(
        xintercept = optimal_cut, linetype = "dashed",
        color = "#D55E00", linewidth = 1.2
      )
    }
    if (criterion == "hazard_ratio") {
      p <- p + ggplot2::geom_hline(yintercept = 1, linetype = "dotted")
    }
    return(p)
  } else if (params$num_cuts == 2) {
    p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$c1, y = .data$c2, fill = .data$stat)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c(name = y_label, option = "viridis") +
      ggplot2::labs(
        title = paste("Systematic Objective Surface:", y_label),
        subtitle = paste("Optimal Matrix Peaks:", paste(round(cutpoint_result$optimal_cuts, 3), collapse = ", ")),
        x = "Cut-point 1 Coordinate", y = "Cut-point 2 Coordinate"
      ) +
      theme_optsurv() +
      ggplot2::theme(panel.grid.major = ggplot2::element_blank())
    return(p)
  } else {
    cli::cli_abort("Surface plotting is restricted to 1 or 2 cuts under systematic grid evaluations.")
  }
}

#' Diagnostic Plot of Schoenfeld Residuals
#'
#' @description
#' High-tier multi-panel diagnostic dashboard tracking the proportional hazards
#' assumption with custom facets per risk cohort stratum.
#'
#' @param x A \code{find_cutpoint} result object.
#' @param ... Unused optional arguments.
#' @return A publication-ready \code{ggplot} canvas frame, or \code{NULL} if the fit fails.
#'
#' @importFrom survival coxph cox.zph Surv
#' @importFrom ggplot2 ggplot aes geom_hline geom_point geom_smooth facet_wrap labs element_blank element_line rel
#' @importFrom cli cli_abort cli_inform
#' @importFrom stats as.formula
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   mock_df <- data.frame(time = 1:30, event = rep(c(0, 1), 15), factor = rnorm(30))
#'   res <- find_cutpoint(
#'      mock_df, "factor", "time", "event",
#'      num_cuts = 1, method = "systematic", quiet = TRUE
#'   )
#'   p <- plot_cutpoint_residuals(res)
#' }
plot_cutpoint_residuals <- function(x, ...) {
  if (is.null(x$optimal_cuts) || any(is.na(x$optimal_cuts))) {
    cli::cli_inform("No valid cut-points mapped; diagnostics skipped.")
    return(invisible(NULL))
  }
  
  data <- x$userdata
  cuts <- sort(x$optimal_cuts)
  num_cuts <- length(cuts)
  
  if (exists("cpp_get_group_assignments", mode = "function") && isTRUE(x$parameters$use_cpp)) {
    grp_vector <- cpp_get_group_assignments(data$factor, cuts)
    data$group <- structure(grp_vector, levels = as.character(1:(num_cuts + 1)), class = "factor")
  } else {
    data$group <- as.factor(findInterval(data$factor, cuts, left.open = TRUE) + 1L)
  }
  
  formula_str <- "survival::Surv(time, event) ~ group"
  covariates <- x$parameters$covariates
  if (!is.null(covariates)) {
    formula_str <- paste(formula_str, "+", paste(covariates, collapse = " + "))
  }
  
  fit <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = data), error = function(e) NULL)
  if (is.null(fit)) {
    cli::cli_inform("Cox fit model failed during residuals generation.")
    return(invisible(NULL))
  }
  
  fit$call$formula <- stats::as.formula(formula_str)
  zph <- tryCatch(survival::cox.zph(fit), error = function(e) NULL)
  
  if (is.null(zph)) {
    cli::cli_abort("Proportional hazards evaluation failed due to a singular model matrix.")
  }
  
  time_vec <- zph$x
  residual_matrix <- zph$y
  col_names <- colnames(residual_matrix)
  
  target_cols = grep("^group", col_names, value = TRUE)
  if (length(target_cols) == 0) {
    cli::cli_inform("No stratified group metrics available for diagnostic modeling.")
    return(invisible(NULL))
  }
  
  long_list <- lapply(target_cols, function(col) {
    data.frame(
      Time = time_vec,
      Residual = residual_matrix[, col],
      Cohort = gsub("group", "Cohort G", col, fixed = TRUE)
    )
  })
  plot_df <- do.call(rbind, long_list)
  
  global_p <- round(zph$table[nrow(zph$table), "p"], 4)
  
  g <- ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$Time, y = .data$Residual)) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.8) +
    ggplot2::geom_point(color = "#2c3e50", alpha = 0.40, size = 1.8, shape = 16) +
    ggplot2::geom_smooth(
      color = "#D55E00", fill = "#D55E00", alpha = 0.15,
      linewidth = 1.2, se = TRUE, method = "loess"
    ) +
    ggplot2::facet_wrap(~ .data$Cohort, scales = "free_y") +
    ggplot2::labs(
      title = "Schoenfeld Residual Diagnostics Dashboard",
      subtitle = paste0("Evaluating Proportional Hazards Assumption | Global Model Fit p = ", global_p),
      x = "Timeline (Follow-up Period Days/Months)",
      y = "Scaled Schoenfeld Residual Metrics"
    ) +
    theme_optsurv()
  
  return(g)
}

#' Plot Landmark Stratification Curves
#'
#' @description Computes and visualizes conditional survival probabilities for patients
#' who survive up to a specified landmark milestone.
#'
#' @param x A \code{find_cutpoint} result object.
#' @param landmark Numeric value indicating the landmark milestone. If NULL,
#'      defaults to 20\% of maximum follow-up data timelines.
#' @param ... Additional arguments passed down to downstream rendering pipelines.
#' @return A \code{ggplot} template or \code{survminer} survival curve asset mapping layout.
#' @export
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   mock_df <- data.frame(
#'      time   = runif(25, 5, 50),
#'      event  = sample(c(0, 1), 25, replace = TRUE),
#'      factor = rnorm(25, 10, 2)
#'   )
#'   res <- find_cutpoint(
#'      mock_df, "factor", "time", "event",
#'      num_cuts = 1, method = "systematic", quiet = TRUE, nmin = 3
#'   )
#'   p <- plot_landmark_stratification(res, landmark = 10)
#' }
plot_landmark_stratification <- function(x, landmark = NULL, ...) {
  df <- x$userdata
  cuts <- sort(x$optimal_cuts)
  num_cuts <- length(cuts)
  
  if (exists("cpp_get_group_assignments", mode = "function") && isTRUE(x$parameters$use_cpp)) {
    grp_vector <- cpp_get_group_assignments(df$factor, cuts)
    df$group <- structure(grp_vector, levels = as.character(1:(num_cuts + 1)), class = "factor")
  } else {
    df$group <- as.factor(findInterval(df$factor, cuts, left.open = TRUE) + 1L)
  }
  
  if (is.null(landmark)) {
    landmark <- round(max(df$time, na.rm = TRUE) * 0.20, 1)
  }
  
  landmark_df <- df[df$time > landmark, ]
  
  if (nrow(landmark_df) < 15) {
    cli::cli_abort("Insufficient sample size remaining (N < 15) to calculate a stable landmark model at time {landmark}.")
  }
  
  landmark_df$time <- landmark_df$time - landmark
  
  cli::cli_alert_info("Generating Landmark Survival Curve for survivors remaining at time milestone: {landmark}")
  
  title_text <- paste("Landmark Survival Analysis (T =", landmark, ")")
  p <- .plot_km_curve(x, landmark_df,
                      title = title_text,
                      xlab = paste("Time Passed Post-Landmark Milestone (Marker Zero =", landmark, ")"), ...
  )
  return(p)
}

#' Plot Cut-point Optimization Stability Surface and S3 Class Interfaces
#'
#' @description
#' Generates a premium, continuous 2D contour surface density topology map tracking the
#' statistical stability of paired discovered cut-points across bootstrap resampling
#' profiles. Automatically handles 1, 2, or multi-cut architectures dynamically.
#'
#' @param validation_result A validation object generated by OptSurvCutR containing
#'     the validation log history dataset.
#' @param main Main title of the chart canvas. Defaults to "Resampling Convergence & Stability Landscape".
#' @param focus_cuts A numeric vector of length 2 specifying which two cut-points to map
#'     if the model contains more than 2 cuts. Defaults to \code{c(1, 2)}.
#' @param x A validation result object for generic dispatch.
#' @param y Unused parameter.
#' @param object A validation result object for summary tracking.
#' @param ... Unused optional dots or arguments passed to printing pipelines.
#' @return A publication-ready \code{ggplot} canvas object displaying stability density bounds.
#'
#' @importFrom ggplot2 ggplot aes geom_density_2d_filled geom_density_2d geom_point labs scale_fill_viridis_d theme element_blank element_line rel geom_density geom_vline
#' @importFrom rlang .data
#' @importFrom cli cli_abort cli_inform cli_h1 cli_text
#' @export
#' @examples
#' mock_val <- list(
#'    bootstrap_distribution = data.frame(Cut_point_1 = rnorm(30, 10, 1)),
#'    original_cuts = 10.2,
#'    parameters = list(predictor = "Biomarker", num_replicates = 30, successful_reps = 30)
#' )
#' p <- plot_validation(mock_val)
plot_validation <- function(validation_result,
                            main = "Resampling Convergence & Stability Landscape",
                            focus_cuts = c(1, 2), ...) {
  log_df <- validation_result$bootstrap_distribution
  
  if (is.null(log_df) || !is.data.frame(log_df)) {
    cli::cli_abort(c(
      "Input must contain a valid statistical log array dataframe inside {.code bootstrap_distribution}.",
      "x" = "Please verify that the validation object was successfully compiled by the bootstrap loop."
    ))
  }
  
  # Extract true original optimal cut-point coordinates
  orig_cuts <- validation_result$original_cuts
  if (is.null(orig_cuts)) orig_cuts <- validation_result$optimal_cuts # Fallback anchor
  
  # Identify available cut-point columns in the dataset
  cut_cols <- grep("^Cut_point_", colnames(log_df), value = TRUE)
  num_available_cuts <- length(cut_cols)
  
  if (num_available_cuts == 0) {
    cli::cli_abort("No columns matching 'Cut_point_X' found in the bootstrap distribution dataset.")
  }
  
  # ===================================================================
  # CONDITIONAL RENDERING MATRIX BASED ON USER CUT ARCHITECTURE
  # ===================================================================
  
  if (num_available_cuts == 1) {
    # --- 1D FALLBACK: Simple single cut-point density ---
    p <- ggplot2::ggplot(log_df, ggplot2::aes(x = .data$Cut_point_1)) +
      ggplot2::geom_density(fill = "#0072B2", alpha = 0.15, color = "#0072B2", linewidth = 1) +
      ggplot2::geom_vline(xintercept = orig_cuts[1], linetype = "dashed", color = "#D55E00", linewidth = 1) +
      ggplot2::labs(
        title = main,
        subtitle = paste0("Optimal Discovered Threshold: ", round(orig_cuts[1], 2)),
        x = paste0(validation_result$parameters$predictor, " Cut-point Coordinate Location"),
        y = "Bootstrap Resampling Density Profile"
      )
  } else {
    # --- 2D CONTOUR ENGINE: For 2, 3, or more cut-points ---
    
    # Safety checks on user input coordinates for multi-cut focus mapping
    if (length(focus_cuts) != 2 || !is.numeric(focus_cuts)) {
      cli::cli_abort("{.arg focus_cuts} must be a numeric vector of length 2 (e.g., c(1, 2) or c(2, 3)).")
    }
    if (any(focus_cuts > num_available_cuts) || any(focus_cuts < 1)) {
      cli::cli_abort("Requested focus cuts [c({focus_cuts[1]}, {focus_cuts[2]})] exceed available model dimensions (Total cuts found: {num_available_cuts}).")
    }
    
    col_x_name <- paste0("Cut_point_", focus_cuts[1])
    col_y_name <- paste0("Cut_point_", focus_cuts[2])
    
    # Extract the paired baseline points for our orange target diamond
    dot_x <- orig_cuts[focus_cuts[1]]
    dot_y <- orig_cuts[focus_cuts[2]]
    
    p <- ggplot2::ggplot(log_df, ggplot2::aes(x = .data[[col_x_name]], y = .data[[col_y_name]])) +
      # Layer 1: Smooth 2D filled density contour ribbons
      ggplot2::geom_density_2d_filled(bins = 10, alpha = 0.95) +
      # Layer 2: White concentric elevation lines
      ggplot2::geom_density_2d(color = "white", alpha = 0.3, linewidth = 0.4, bins = 10) +
      # Layer 3: Target Baseline Crosshair Anchor (Tidy-evaluation compliant)
      ggplot2::geom_point(
        data = data.frame(x = dot_x, y = dot_y),
        ggplot2::aes(x = .data$x, y = .data$y), inherit.aes = FALSE,
        shape = 23, size = 4.5, fill = "#D55E00", color = "white", stroke = 1.2
      ) +
      ggplot2::scale_fill_viridis_d(name = "Discovery Density", option = "mako") +
      ggplot2::labs(
        title = main,
        subtitle = paste0(
          "Evaluating Cut ", focus_cuts[1], " vs. Cut ", focus_cuts[2],
          " | Original Coordinates: [", round(dot_x, 2), ", ", round(dot_y, 2), "]"
        ),
        x = paste0(validation_result$parameters$predictor, " Cut-point ", focus_cuts[1], " Threshold"),
        y = paste0(validation_result$parameters$predictor, " Cut-point ", focus_cuts[2], " Threshold")
      )
  }
  
  # Apply uniform package aesthetic parameters
  p <- p +
    theme_optsurv() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "grey70", linewidth = 0.6),
      legend.position = "right",
      legend.title = ggplot2::element_text(size = ggplot2::rel(0.85), face = "bold"),
      legend.text = ggplot2::element_text(size = ggplot2::rel(0.80))
    )
  
  return(p)
}

# --- S3 OVERRIDE INTERFACES REGULATING CLASSES ---

#' @rdname plot_validation
#' @export
plot.validate_cutpoint_result <- function(x, y, ...) {
  plot_validation(validation_result = x, ...)
}

#' @rdname plot_validation
#' @export
print.validate_cutpoint_result <- function(x, ...) {
  cli::cli_h1("Validation Results Summary")
  cli::cli_text("Replicates run: {.val {x$parameters$num_replicates}} ({.val {x$parameters$successful_reps}} successful)")
  cli::cli_inform("Use summary() to view confidence intervals matrix profiles.")
  invisible(x)
}

#' @rdname plot_validation
#' @export
summary.validate_cutpoint_result <- function(object, ...) {
  cli::cli_h1("Validation Stability Analytics")
  print(object$confidence_intervals)
  invisible(object)
}

# --- BACKGROUND INTERIOR ENGINE PLOTTING MODULES ---

.plot_km_curve <- function(x, df, title = "Kaplan-Meier Survival Estimation",
                           xlab = "Follow-up Time", ylab = "Overall Survival Probability", ...) {
  if (!requireNamespace("survminer", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg survminer} is required to render outcome tracking charts.")
  }
  
  target_formula <- stats::as.formula("survival::Surv(time, event) ~ group")
  fit_km <- survival::survfit(target_formula, data = df)
  fit_km$call$formula <- target_formula
  
  # Count the actual number of categorical risk strata present in the model fit
  num_strata <- length(names(fit_km$strata))
  
  # Generate a publication-ready color spectrum that scales dynamically to match num_strata
  dynamic_palette <- if (num_strata <= 4) {
    c("#0072B2", "#D55E00", "#009E73", "#CC79A7")[1:num_strata]
  } else {
    grDevices::colorRampPalette(c("#0072B2", "#009E73", "#D55E00", "#CC79A7"))(num_strata)
  }
  
  p <- survminer::ggsurvplot(fit_km,
                             data = df, title = title, xlab = xlab, ylab = ylab,
                             palette = dynamic_palette,
                             pval = TRUE, ggtheme = theme_optsurv(), ...
  )
  
  if (!is.null(p$plot)) p$plot <- p$plot + theme_optsurv()
  if (!is.null(p$table)) p$table <- p$table + theme_optsurv()
  if (!is.null(p$ncensor.plot)) p$ncensor.plot <- p$ncensor.plot + theme_optsurv()
  
  return(p)
}

.plot_density_cuts <- function(x, df, title = "Predictor Cohort Distribution Map", ...) {
  p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$factor)) +
    ggplot2::geom_density(fill = "#0072B2", alpha = 0.15, color = "#0072B2", linewidth = 1) +
    ggplot2::labs(
      title = title, subtitle = paste("Marker:", x$parameters$predictor),
      x = paste(x$parameters$predictor, "(Continuous Range)"), y = "Population Density Profile"
    ) +
    theme_optsurv()
  
  stagger_vjust <- c(1.5, 3.8, 6.1, 8.4)
  idx <- 1
  
  for (cut_val in x$optimal_cuts) {
    current_vjust <- stagger_vjust[((idx - 1) %% length(stagger_vjust)) + 1]
    
    p <- p +
      ggplot2::geom_vline(
        xintercept = cut_val, linetype = "dashed",
        color = "#D55E00", linewidth = 0.8, alpha = 0.5
      ) +
      ggplot2::annotate("text",
                        x = cut_val, y = Inf,
                        label = paste0("Cut ", idx, ": ", round(cut_val, 2)),
                        vjust = current_vjust, hjust = -0.1,
                        color = "#2c3e50", fontface = "bold", size = 3.5
      )
    idx <- idx + 1
  }
  return(p)
}

#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Displays secondary effect sizes alongside corresponding model covariates.
#' @importFrom ggplot2 ggplot aes geom_blank theme geom_vline geom_errorbar geom_point geom_col scale_fill_manual labs element_blank margin element_line element_rect geom_text scale_color_manual xlim rel
#' @importFrom patchwork wrap_plots plot_layout
#' @importFrom stats as.formula aggregate relevel symnum
#' @importFrom rlang .data
.plot_hr_forest <- function(x, df, reference_group = "G1", main = "Adjusted Clinical Risk Profile", ...) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg survival} is required to compute forest metrics.")
  }
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    cli::cli_abort("Package {.pkg patchwork} is required for space-optimized ribbon forest plots.")
  }
  
  covariates <- x$parameters$covariates
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""
  
  # Reconstruct biomarker risk cohorts from optimal cuts
  num_cuts <- length(x$optimal_cuts)
  df$group <- cut(df$factor, breaks = c(-Inf, sort(x$optimal_cuts), Inf), labels = paste0("G", 1:(num_cuts + 1)))
  
  if (reference_group %in=% levels(df$group)) {
    df$group <- stats::relevel(df$group, ref = reference_group)
  }
  
  # 1. Fit Cox Regression Model
  formula_str <- paste("survival::Surv(time, event) ~ group", cov_part)
  fit <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = df), error = function(e) NULL)
  if (is.null(fit)) cli::cli_abort("Could not compile forest plots: Cox regression alignment failure.")
  
  sum_cox <- summary(fit)
  coefs <- as.data.frame(sum_cox$coefficients)
  conf <- as.data.frame(sum_cox$conf.int)
  
  # Identify model coefficients
  biomarker_rows <- grep("^group", rownames(coefs))
  covariate_rows <- which(!seq_len(nrow(coefs)) %in% biomarker_rows)
  
  if (length(biomarker_rows) == 0) cli::cli_abort("No valid categorical group contrasts found to chart.")
  
  # 2. Calculate Biomarker Cohort Sample Sizes (N & Events)
  counts_df <- as.data.frame(table(df$group))
  names(counts_df) <- c("g_id", "N")
  event_df <- stats::aggregate(event ~ group, data = df, sum)
  names(event_df) <- c("g_id", "Events")
  km_sum <- merge(counts_df, event_df, by = "g_id")
  rownames(km_sum) <- km_sum$g_id
  
  # 3. Compile Master Dataset for the Left Panel (Forest Rows)
  forest_list <- list()
  
  # A. Inject Reference Group Row Explicitly
  ref_label <- paste0("Cohort ", reference_group, " (Reference)")
  forest_list[[1]] <- data.frame(
    Variable  = ref_label,
    Type      = "Biomarker",
    HR        = 1.0,
    Lower     = 1.0,
    Upper     = 1.0,
    HRText    = "Reference",
    SortOrder = 1,
    g_id      = reference_group
  )
  
  # B. Inject Stratified Biomarker Rows
  idx <- 2
  for (r_idx in biomarker_rows) {
    rname <- rownames(coefs)[r_idx]
    g_id <- gsub("group", "", rname, fixed = TRUE)
    pval <- coefs[rname, "Pr(>|z|)"]
    p_text <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))
    
    hr_val <- conf[r_idx, "exp(coef)"]
    low_val <- conf[r_idx, "lower .95"]
    up_val <- conf[r_idx, "upper .95"]
    
    forest_list[[idx]] <- data.frame(
      Variable  = paste0("Cohort ", g_id, " (", p_text, ")"),
      Type      = "Biomarker",
      HR        = hr_val,
      Lower     = low_val,
      Upper     = up_val,
      HRText    = paste0("HR = ", round(hr_val, 2), " (", round(low_val, 2), "-", round(up_val, 2), ")"),
      SortOrder = idx,
      g_id      = g_id
    )
    idx <- idx + 1
  }
  
  # C. Inject Clinical Adjusters Covariate Rows (Dynamic Coloring Category)
  if (length(covariate_rows) > 0) {
    for (r_idx in covariate_rows) {
      rname <- rownames(coefs)[r_idx]
      pval <- coefs[r_idx, "Pr(>|z|)"]
      p_text <- if (pval < 0.001) "p < 0.001" else paste0("p = ", round(pval, 3))
      
      hr_val <- conf[r_idx, "exp(coef)"]
      low_val <- conf[r_idx, "lower .95"]
      up_val <- conf[r_idx, "upper .95"]
      
      forest_list[[idx]] <- data.frame(
        Variable  = paste0(rname, " (", p_text, ")"),
        Type      = "Covariate",
        HR        = hr_val,
        Lower     = low_val,
        Upper     = up_val,
        HRText    = paste0("HR = ", round(hr_val, 2), " (", round(low_val, 2), "-", round(up_val, 2), ")"),
        SortOrder = idx,
        g_id      = NA
      )
      idx <- idx + 1
    }
  }
  
  forest_df <- do.call(rbind, forest_list)
  forest_df$Variable <- factor(forest_df$Variable, levels = rev(forest_df$Variable[order(forest_df$SortOrder)]))
  
  # 4. Assemble Right Panel Dataset (Guaranteed level-safe binding match)
  bar_list <- list()
  b_idx <- 1
  
  for (f_row in seq_len(nrow(forest_df))) {
    current_gid <- forest_df$g_id[f_row]
    current_var <- forest_df$Variable[f_row]
    
    if (!is.na(current_gid)) {
      bar_list[[b_idx]] <- data.frame(
        Variable = rep(current_var, 2),
        Value    = c(km_sum[current_gid, "N"] - km_sum[current_gid, "Events"], km_sum[current_gid, "Events"]),
        Metric   = c("Censored", "Events")
      )
    } else {
      bar_list[[b_idx]] <- data.frame(
        Variable = rep(current_var, 2),
        Value    = c(0, 0),
        Metric   = c("Censored", "Events")
      )
    }
    b_idx <- b_idx + 1
  }
  
  bar_df <- do.call(rbind, bar_list)
  bar_df$Variable <- factor(bar_df$Variable, levels = levels(forest_df$Variable))
  bar_df$Metric <- factor(bar_df$Metric, levels = c("Censored", "Events"))
  
  max_upper <- max(forest_df$Upper, na.rm = TRUE)
  x_limit_upper <- max_upper * 1.05
  
  # ===================================================================
  # PANE 1: THE REFACTORED FOREST DATA CANVAS
  # ===================================================================
  p_forest <- ggplot2::ggplot(forest_df, ggplot2::aes(y = .data$Variable)) +
    theme_optsurv() +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor.y = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_line(color = "grey88", linewidth = 0.5),
      panel.grid.minor.x = ggplot2::element_line(color = "grey94", linewidth = 0.25),
      legend.position = "none",
      plot.margin = ggplot2::margin(5, 0, 5, 5)
    ) +
    ggplot2::geom_vline(xintercept = 1.0, linetype = "dotted", linewidth = 0.9, color = "gray30") +
    # --- Modern ggplot2 Horizontal Layout ---
    ggplot2::geom_errorbar(
      ggplot2::aes(xmin = .data$Lower, xmax = .data$Upper, color = .data$Type),
      width = 0.12,
      linewidth = 0.7,
      orientation = "y"
    ) +
    ggplot2::geom_point(ggplot2::aes(x = .data$HR, color = .data$Type), shape = 16, size = 3.5) +
    ggplot2::scale_color_manual(values = c("Biomarker" = "black", "Covariate" = "#7f8c8d")) +
    ggplot2::xlim(0, x_limit_upper) +
    ggplot2::labs(title = main, subtitle = NULL, x = "Hazard Ratio (95% CI)", y = "")
  
  # ===================================================================
  # PANE 2: THE VARIABLE-SCALE SIDEBAR RIBBON (1:10 Ratio Scale)
  # ===================================================================
  p_ribbon <- ggplot2::ggplot(bar_df, ggplot2::aes(y = .data$Variable, x = .data$Value, fill = .data$Metric)) +
    theme_optsurv() +
    ggplot2::theme(
      panel.grid.major = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.x = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 0),
      legend.position = "bottom"
    ) +
    ggplot2::geom_col(width = 0.35, position = "stack") +
    ggplot2::scale_fill_manual(name = "", values = c("Censored" = "#0072B2", "Events" = "#D55E00"))
  
  combined_layout <- patchwork::wrap_plots(p_forest, p_ribbon, ncol = 2, widths = c(10, 1)) +
    patchwork::plot_layout(guides = "collect") &
    ggplot2::theme(legend.position = "bottom")
  
  return(combined_layout)
}

.plot_genetic_trajectory <- function(x, ...) {
  cli::cli_inform("Trajectory tracking plots (`type = 'trajectory'`) are reserved for downstream optimization tracking structures in evolutionary models.")
  return(invisible(NULL))
}
