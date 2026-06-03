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
  expect_true(inherits(t$panel.grid.minor, "element_blank"))
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
  expect_equal(p_hr$labels$y, "Hazard Ratio (HR)")

  res_pv <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "p_value", nmin = 1, quiet = TRUE)
  ))
  skip_if(is.null(res_pv) || all(is.na(res_pv$optimal_cuts)), "Systematic p-value search failed.")

  p_pv <- plot_optimisation_curve(res_pv)
  expect_s3_class(p_pv, "ggplot")
  expect_equal(p_pv$labels$y, "P-value")
})

#' @srrstats {G2.0} Tests input validation for plotting function.
test_that("plot_optimisation_curve throws errors for invalid input", {
  res_valid <- valid_fc_result_for_boot

  expect_error(plot_optimisation_curve(list()), regexp = "Input must be an object from the")

  res_genetic <- res_valid
  res_genetic$parameters$method <- "genetic"
  expect_error(plot_optimisation_curve(res_genetic), regexp = "only for `method = \"systematic\"`")

  res_2_cuts <- res_valid
  res_2_cuts$parameters$num_cuts <- 2
  expect_error(plot_optimisation_curve(res_2_cuts), regexp = "only for `num_cuts = 1`")

  res_no_stats <- res_valid
  res_no_stats$all_stats <- NULL
  expect_error(plot_optimisation_curve(res_no_stats), regexp = "must have a valid `all_stats`")

  res_valid$all_stats <- data.frame()
  expect_error(plot_optimisation_curve(res_valid), regexp = "`all_stats` is empty")
})

#' @srrstats {RE6.3} Tests diagnostic Schoenfeld residual plot.
test_that("plot_cutpoint_residuals produces Schoenfeld residual plot (RE6.3)", {
  skip_if(any(is.na(valid_fc_result_for_boot$optimal_cuts)))

  p <- plot_cutpoint_residuals(valid_fc_result_for_boot)
  expect_type(p, "list")
  expect_true(length(p) > 0)
  expect_true(all(vapply(p, inherits, "ggplot", FUN.VALUE = logical(1))))

  titles <- vapply(p, function(x) x$labels$title %||% "", character(1))
  expect_true(any(grepl("Schoenfeld", titles, ignore.case = TRUE)))

  last_subtitle <- p[[length(p)]]$labels$subtitle %||% ""
  expect_match(last_subtitle, "Global", all = FALSE)
})

test_that("Coverage: plot_cutpoint_residuals handles data failure modes", {
  # Pathological events (all 0) - coxph survives, but cox.zph crashes
  res_no_event <- valid_fc_result_for_boot
  res_no_event$userdata$event <- 0
  expect_error(
    plot_cutpoint_residuals(res_no_event),
    regexp = "Proportional hazards test failed"
  )

  # Corrupted data types - coxph crashes immediately
  res_corrupt <- valid_fc_result_for_boot
  res_corrupt$userdata$event <- as.character(res_corrupt$userdata$event)
  expect_message(
    plot_cutpoint_residuals(res_corrupt),
    regexp = "Cox model failed"
  )
})

test_that("plot_time_dependent_auc successfully generates an AUC plot", {
  skip_if_not_installed("timeROC")
  skip_if(any(is.na(valid_fc_result_for_boot$optimal_cuts)))

  p_auc <- plot_time_dependent_auc(valid_fc_result_for_boot)
  expect_s3_class(p_auc, "ggplot")
  expect_equal(p_auc$labels$title, "Time-Dependent AUC")
  expect_equal(p_auc$labels$y, "Area Under the Curve (AUC)")
})

test_that("plot_time_dependent_auc handles pathological data (insufficient events)", {
  skip_if_not_installed("timeROC")

  res_few_events <- valid_fc_result_for_boot
  res_few_events$userdata$event <- 0
  res_few_events$userdata$event[1:5] <- 1

  expect_error(plot_time_dependent_auc(res_few_events), regexp = "Not enough events for stable time-dependent AUC")
})

# --- 3. S3 Router (plot.find_cutpoint) ---

#' @srrstats {RE6.0} Tests S3 plot method.
test_that("plot.find_cutpoint generates all plot types and routes correctly", {
  skip_if_not_installed("broom")
  skip_if(any(is.na(valid_fc_result_for_boot$optimal_cuts)))

  res_valid <- valid_fc_result_for_boot

  # Outcome Plot
  p_outcome <- suppressWarnings(plot(res_valid, type = "outcome"))
  expect_s3_class(p_outcome$plot, "ggplot")

  # Distribution Plot
  p_dist <- plot(res_valid, type = "distribution")
  expect_s3_class(p_dist, "ggplot")
  expect_equal(p_dist$labels$y, "Density")

  # Forest Plot
  p_forest <- plot(res_valid, type = "forest")
  expect_s3_class(p_forest, "ggplot")
  expect_equal(p_forest$labels$x, "Hazard Ratio (95% CI)")

  # Diagnostic Plot (Routing)
  p_diag <- suppressWarnings(plot(res_valid, type = "diagnostic"))
  expect_true(inherits(p_diag, "ggcoxzph") || is.list(p_diag))

  # AUC Plot (Routing)
  if (requireNamespace("timeROC", quietly = TRUE)) {
    p_auc <- plot(res_valid, type = "auc")
    expect_s3_class(p_auc, "ggplot")
  }

  # Dashboard Plot (Routing)
  if (requireNamespace("patchwork", quietly = TRUE)) {
    p_dash <- suppressWarnings(plot(res_valid, type = "all"))
    expect_s3_class(p_dash, "patchwork")
  }

  # Escape Hatch (Data extraction)
  raw_df <- plot(res_valid, return_data = TRUE)
  expect_true(is.data.frame(raw_df))
  expect_true("group" %in% colnames(raw_df))
})

test_that("plot.find_cutpoint references and handles Cox model failure", {
  skip_if_not_installed("broom")

  # 2-Cut Reference Group Test
  res_2_cuts <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event", num_cuts = 2, method = "genetic", max.generations = 5, quiet = TRUE)
  ))
  skip_if(any(is.na(res_2_cuts$optimal_cuts)))

  p_forest_g2 <- plot(res_2_cuts, type = "forest", reference_group = "G2")
  expect_s3_class(p_forest_g2, "ggplot")
  expect_match(p_forest_g2$labels$subtitle, "G2")

  # Invalid Reference Group
  expect_message(plot(res_2_cuts, type = "forest", reference_group = "INVALID"), regexp = "Defaulting to")

  # S3 Cox Model Failure
  local_mocked_bindings("coxph" = function(...) NULL, .package = "survival")
  expect_message(plot(valid_fc_result_for_boot, type = "forest"), regexp = "Could not fit Cox model for forest plot")
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
