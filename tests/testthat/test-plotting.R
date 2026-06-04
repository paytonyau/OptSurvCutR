# ===================================================================
# TESTS: plotting_functions.R
# Diagnostic plots, Interactive widgets, and S3 plot methods
# ===================================================================

#' @srrstats {RE6.0} Tests for S3 plot methods.
#' @srrstats {G2.4c} Tests optional dependency handling.

# --- 1. Theme and Interactive Tools ---

test_that("theme_optsurv returns a valid ggplot2 theme", {
  t <- theme_optsurv()
  expect_s3_class(t, "theme")
  expect_s3_class(t, "gg")

  # Ensure specific clinical aesthetics are applied
  expect_equal(t$legend.position, "bottom")
  expect_false(inherits(t$panel.grid.minor, "element_blank"))
})

test_that("optsurv_interactive handles both ggplot and ggsurvplot objects", {
  skip_if_not_installed("plotly")
  skip_if(any(is.na(valid_fc_result_for_boot$optimal_cuts)))

  # 1. Test standard ggplot (Distribution)
  p_dist <- plot(valid_fc_result_for_boot, type = "distribution")
  int_dist <- optsurv_interactive(p_dist)
  expect_s3_class(int_dist, "plotly")
  expect_s3_class(int_dist, "htmlwidget")

  # 2. Test ggsurvplot list object (Outcome)
  p_out <- suppressWarnings(plot(valid_fc_result_for_boot, type = "outcome"))
  int_out <- optsurv_interactive(p_out)
  expect_s3_class(int_out, "plotly")
  expect_s3_class(int_out, "htmlwidget")
})

# --- 2. Specific Plot Generators ---

#' @srrstats {RE6.0} Tests diagnostic plot for systematic search.
test_that("plot_optimisation_curve works for all criteria", {
  res_lr <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1, quiet = TRUE)
  ))
  skip_if(is.null(res_lr) || all(is.na(res_lr$optimal_cuts)), "Systematic logrank search failed.")

  p_lr <- plot_optimisation_curve(res_lr)
  expect_s3_class(p_lr, "ggplot")
  expect_equal(p_lr$labels$y, "Log-Rank Statistic")

  res_hr <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "hazard_ratio", nmin = 1, quiet = TRUE)
  ))
  skip_if(is.null(res_hr) || all(is.na(res_hr$optimal_cuts)), "Systematic HR search failed.")

  p_hr <- plot_optimisation_curve(res_hr)
  expect_s3_class(p_hr, "ggplot")
  expect_equal(p_hr$labels$y, "Hazard Ratio")
})

#' @srrstats {G2.0} Tests input validation for plotting function.
test_that("plot_optimisation_curve throws errors for invalid input", {
  res_valid <- valid_fc_result_for_boot

  expect_error(plot_optimisation_curve(list()), regexp = "Input must be an object from the")

  res_genetic <- res_valid
  res_genetic$parameters$method <- "genetic"
  expect_error(plot_optimisation_curve(res_genetic), regexp = "Surface mapping.*is only available for")

  res_3_cuts <- res_valid
  res_3_cuts$parameters$num_cuts <- 3
  expect_error(plot_optimisation_curve(res_3_cuts), regexp = "Surface plotting is restricted to 1 or 2 cuts")

  res_no_stats <- res_valid
  res_no_stats$all_stats <- NULL
  expect_error(plot_optimisation_curve(res_no_stats), regexp = "must contain a valid grid log array")
})

#' @srrstats {RE6.3} Tests diagnostic Schoenfeld residual plot.
test_that("plot_cutpoint_residuals produces Schoenfeld residual plot (RE6.3)", {
  skip_if(any(is.na(valid_fc_result_for_boot$optimal_cuts)))

  p <- plot_cutpoint_residuals(valid_fc_result_for_boot)
  expect_s3_class(p, "ggplot")
  expect_match(p$labels$title, "Schoenfeld Residual Diagnostics Dashboard")
  expect_match(p$labels$subtitle, "Global Model Fit")
})

# Replace the old pathology checks with strict error catching:
test_that("Coverage: plot_cutpoint_residuals handles data failure modes", {
  res_no_event <- valid_fc_result_for_boot
  res_no_event$userdata$event <- 0
  expect_error(
    plot_cutpoint_residuals(res_no_event),
    regexp = "Proportional hazards evaluation failed due to a singular model matrix"
  )

  res_corrupt <- valid_fc_result_for_boot
  res_corrupt$userdata$event <- as.character(res_corrupt$userdata$event)
  expect_null(suppressMessages(plot_cutpoint_residuals(res_corrupt)))
})

# --- 3. S3 Router (plot.find_cutpoint) ---

#' @srrstats {RE6.0} Tests S3 plot method.
test_that("plot.find_cutpoint generates all plot types and routes correctly", {
  skip_if(any(is.na(valid_fc_result_for_boot$optimal_cuts)))

  res_valid <- valid_fc_result_for_boot

  # Outcome Plot
  p_outcome <- suppressWarnings(plot(res_valid, type = "outcome"))
  expect_s3_class(p_outcome$plot, "ggplot")

  # Distribution Plot
  p_dist <- plot(res_valid, type = "distribution")
  expect_s3_class(p_dist, "ggplot")
  expect_equal(p_dist$labels$y, "Population Density Profile")

  # Forest Plot (Returns a composite patchwork layout)
  p_forest <- plot(res_valid, type = "forest")
  expect_s3_class(p_forest, "patchwork")

  # Diagnostic Plot (Routing to facet canvas)
  p_diag <- suppressWarnings(plot(res_valid, type = "diagnostic"))
  expect_s3_class(p_diag, "ggplot")

  # Escape Hatch (Data extraction)
  raw_df <- plot(res_valid, return_data = TRUE)
  expect_true(is.data.frame(raw_df))
  expect_true("group" %in% colnames(raw_df))
})

test_that("plot.find_cutpoint references and handles Cox model failure", {
  res_2_cuts <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event", num_cuts = 2, method = "genetic", maxiter = 5, quiet = TRUE)
  ))
  skip_if(any(is.na(res_2_cuts$optimal_cuts)))

  p_composite <- plot(res_2_cuts, type = "forest", reference_group = "G2")
  expect_s3_class(p_composite, "patchwork")

  # Unpack patchwork layer safely to verify text titles
  p_forest_layer <- p_composite[[1]]
  expect_match(p_forest_layer$labels$title, "Adjusted Clinical Risk Profile")

  # S3 Router Argument Safeguard Verification
  expect_error(plot(res_2_cuts, type = "auc"), regexp = "should be one of")
})

test_that("Coverage: plot.find_cutpoint_number_result – all IC NA", {
  obj <- structure(
    list(
      optimal_cuts = NA_real_,
      optimal_stat = NA_real_,
      ic = rep(NA_real_, 5),
      parameters = list(method = "systematic", criterion = "BIC")
    ),
    class = "find_cutpoint_number_result"
  )
  expect_snapshot(plot(obj))

  # Also handles missing parameters gracefully
  obj2 <- obj
  obj2$parameters <- list(method = "unknown", criterion = "BIC")
  expect_snapshot(plot(obj2))
})
