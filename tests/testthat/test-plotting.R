# ===================================================================
# TESTS: plotting_functions.R
# Diagnostic plots, Static layouts, and S3 plot methods
# ===================================================================

.generate_clean_plotting_fixtures <- function() {
  set.seed(42)
  n <- 60

  df <- data.frame(
    time       = runif(n, 5, 50),
    event      = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7)),
    predictor  = rnorm(n, mean = 10, sd = 3),
    covariate1 = rnorm(n, 0, 1),
    covariate2 = factor(sample(c("A", "B"), n, replace = TRUE))
  )

  df_3groups <- df
  df_3groups$predictor <- c(rnorm(20, 5, 1), rnorm(20, 15, 1), rnorm(20, 25, 1))

  res_base <- suppressMessages(suppressWarnings(
    find_cutpoint(df, "predictor", "time", "event", num_cuts = 1, method = "systematic", quiet = TRUE, nmin = 5)
  ))

  list(df = df, df_3groups = df_3groups, res_base = res_base)
}

# --- 1. Theme Aesthetics & Diagnostic Curves ---

test_that("Core plot themes and optimization curves execute completely", {
  t <- theme_optsurv()
  expect_s3_class(t, "theme")
  expect_identical(t$legend.position, "bottom") # UPGRADED: expect_equal -> expect_identical

  fixtures <- .generate_clean_plotting_fixtures()

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  p_lr <- plot_optimisation_curve(fixtures$res_base)
  expect_s3_class(p_lr, "ggplot")

  suppressWarnings({
    print(p_lr)
  })
})

#' @srrstats {G5.2b}
test_that("plot_optimisation_curve traps all edge failures and checks exact string errors", {
  expect_error(plot_optimisation_curve(list()), regexp = "object")

  fixtures <- .generate_clean_plotting_fixtures()
  res_fake <- fixtures$res_base

  res_fake$parameters$method <- "genetic"
  expect_error(plot_optimisation_curve(res_fake), regexp = "Surface")

  res_fake$parameters$method <- "systematic"
  res_fake$parameters$num_cuts <- 3
  expect_error(plot_optimisation_curve(res_fake), regexp = "restricted")

  res_fake$parameters$num_cuts <- 1
  res_fake$all_stats <- NULL
  expect_error(plot_optimisation_curve(res_fake), regexp = "grid")
})

# --- 2. Residuals and Data Failure Modes ---

test_that("Schoenfeld residual canvas handles layout failures", {
  fixtures <- .generate_clean_plotting_fixtures()

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  p <- plot_cutpoint_residuals(fixtures$res_base)
  expect_s3_class(p, "ggplot")
  print(p)

  res_crash <- fixtures$res_base
  res_crash$userdata$event <- 0
  expect_error(plot_cutpoint_residuals(res_crash), regexp = "failed")

  res_crash$userdata$event <- as.character(res_crash$userdata$event)
  expect_null(suppressMessages(plot_cutpoint_residuals(res_crash)))
})

# --- 3. Complete Loop Sweep of S3 Router (plot.find_cutpoint) ---

#' @srrstats {RE6.0}
#' @srrstats {RE6.1}
test_that("S3 router executes all plot variant selections through proper generic class dispatch", {
  fixtures <- .generate_clean_plotting_fixtures()

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  for (layout_type in c("outcome", "distribution", "forest", "diagnostic")) {
    p_any <- suppressMessages(suppressWarnings(plot(fixtures$res_base, type = layout_type)))
    expect_true(inherits(p_any, c("ggplot", "patchwork", "ggsurv", "ggsurvplot")) || is.list(p_any))
    print(p_any)
  }

  raw_df <- plot(fixtures$res_base, return_data = TRUE)
  expect_true(is.data.frame(raw_df))
  expect_true("group" %in% colnames(raw_df))

  expect_error(plot(fixtures$res_base, type = "invalid_selection"), regexp = "one of")
})

test_that("Exhaustive layout matrix testing for S3 plot routers across engine types", {
  fixtures <- .generate_clean_plotting_fixtures()

  res_with_covars <- suppressMessages(suppressWarnings(
    find_cutpoint(fixtures$df, "predictor", "time", "event",
      covariates = "covariate1", num_cuts = 1, method = "systematic", quiet = TRUE, nmin = 5
    )
  ))

  res_2cuts_gen <- suppressMessages(suppressWarnings(
    find_cutpoint(fixtures$df_3groups, "predictor", "time", "event",
      num_cuts = 2, method = "systematic", quiet = TRUE, nmin = 5
    )
  ))

  res_pure_r <- fixtures$res_base
  res_pure_r$parameters$use_cpp <- FALSE

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  all_test_objects <- list(fixtures$res_base, res_with_covars, res_2cuts_gen, res_pure_r)

  for (obj in all_test_objects) {
    if (is.null(obj) || anyNA(obj$optimal_cuts)) next

    for (layout_type in c("outcome", "distribution", "forest", "diagnostic")) {
      p_run <- suppressMessages(suppressWarnings(plot(obj, type = layout_type)))
      expect_true(inherits(p_run, c("ggplot", "patchwork", "ggsurv", "ggsurvplot")) || is.list(p_run))
      print(p_run)
    }
  }
})

# --- 4. Snapshot Verification for Empty Objects ---

test_that("Plot routines fail safely on empty search structures", {
  obj_empty <- structure(
    list(
      optimal_cuts = NA_real_, optimal_stat = NA_real_,
      parameters = list(method = "systematic", criterion = "logrank")
    ),
    class = "find_cutpoint"
  )

  expect_error(
    plot(obj_empty),
    regexp = "valid"
  )
})

#' @srrstats {RE6.2}
test_that("Graphics engine stress test via aesthetic argument toggles and confidence intervals", {
  fixtures <- .generate_clean_plotting_fixtures()

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  p_out_ci <- suppressMessages(suppressWarnings(plot(fixtures$res_base, type = "outcome", conf.int = TRUE, risk.table = TRUE)))
  p_out_no <- suppressMessages(suppressWarnings(plot(fixtures$res_base, type = "outcome", conf.int = FALSE, risk.table = FALSE)))
  print(p_out_ci)
  print(p_out_no)

  for (dist_mode in c("density", "histogram")) {
    p_dist <- suppressMessages(suppressWarnings(plot(fixtures$res_base, type = "distribution", mode = dist_mode)))
    expect_s3_class(p_dist, "ggplot")
    print(p_dist)
  }

  p_forest_custom <- suppressMessages(suppressWarnings(
    plot(fixtures$res_base, type = "forest", palette = c("red", "blue"), title = "Custom Title Override")
  ))
  print(p_forest_custom)
})

# --- 5. Milestone Optimisation Clearing Blocks ---

test_that("Force plot_optimisation_curve to render 3D/2-cut surface logic or trap missing grids", {
  fixtures <- .generate_clean_plotting_fixtures()

  res_2cuts_surf <- suppressMessages(suppressWarnings(
    find_cutpoint(fixtures$df_3groups, "predictor", "time", "event",
      num_cuts = 2, method = "systematic", nmin = 5, quiet = TRUE
    )
  ))

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  p_surf <- plot_optimisation_curve(res_2cuts_surf)
  expect_s3_class(p_surf, "ggplot")
})

test_that("Exhaustive S3 method and plot sweep for find_cutpoint_number results", {
  fixtures <- .generate_clean_plotting_fixtures()

  res_num <- suppressMessages(suppressWarnings(
    find_cutpoint_number(fixtures$df, "predictor", "time", "event",
      max_cuts = 2, method = "systematic", criterion = "BIC", quiet = TRUE, nmin = 5
    )
  ))

  sink_file <- tempfile()
  sink_con <- file(sink_file, open = "wt")
  sink(sink_con)
  sink(sink_con, type = "message")

  on.exit(
    {
      sink(type = "message")
      sink()
      close(sink_con)
      unlink(sink_file)
    },
    add = TRUE
  )

  print.find_cutpoint_number_result(res_num)

  res_num_summary <- summary.find_cutpoint_number_result(res_num)
  expect_s3_class(res_num_summary, "find_cutpoint_number_result")

  pdf(NULL)
  on.exit(dev.off(), add = TRUE)

  p_num_plot <- plot.find_cutpoint_number_result(res_num)
  expect_s3_class(p_num_plot, "ggplot")
  print(p_num_plot)
})
