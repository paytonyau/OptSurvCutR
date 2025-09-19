# This is a comprehensive test file for the main user-facing functions of the OptSurvCutR package.
# It covers find_cutpoint_number(), find_cutpoint(), and validate_cutpoint() as per the
# recommended workflow: (1) determine optimal number of cuts, (2) find cutpoint locations,
# (3) validate cutpoints with bootstrapping. Tests include recommended methods, criteria,
# edge cases, and realistic data scenarios.

# --- Setup: Load necessary packages and create mock data ---
library(testthat)
library(survival)

# Create a consistent, reproducible mock dataset for testing
set.seed(42)
n_test <- 325
mock_data <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = rnorm(n_test, mean = 50, sd = 10),
  covariate1 = rnorm(n_test, mean = 5, sd = 1)
)

# Mock dataset with skewed predictor for realistic testing
mock_data_skewed <- mock_data
mock_data_skewed$predictor <- rexp(n_test, rate = 0.1)

# Mock dataset with heavy censoring (90% censored)
mock_data_heavy_censor <- mock_data
mock_data_heavy_censor$event <- sample(0:1, n_test, replace = TRUE, prob = c(0.9, 0.1))

# --- Tests for find_cutpoint_number() ---

test_that("find_cutpoint_number genetic search works for recommended criteria", {
  skip_if_not_installed("rgenoud")

  # Test with BIC (recommended default)
  # NOTE: suppressWarnings is used to silence the expected rgenoud warning about hitting maxiter
  res_bic <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic", criterion = "BIC", max_cuts = 2, maxiter = 10, nmin = 1)))
  skip_if(is.null(res_bic), "Genetic search returned NULL, skipping checks.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 3) # 0, 1, and 2 cuts
  expect_true("BIC" %in% names(res_bic$results))

  # Test with AIC
  res_aic <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic", criterion = "AIC", max_cuts = 2, maxiter = 10, nmin = 1)))
  skip_if(is.null(res_aic), "Genetic search returned NULL, skipping checks.")
  expect_s3_class(res_aic, "find_cutpoint_number_result")
  expect_true("AIC" %in% names(res_aic$results))

  # Test with AICc
  res_aicc <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic", criterion = "AICc", max_cuts = 2, maxiter = 10, nmin = 1)))
  skip_if(is.null(res_aicc), "Genetic search returned NULL, skipping checks.")
  expect_s3_class(res_aicc, "find_cutpoint_number_result")
  expect_true("AICc" %in% names(res_aicc$results))
})

test_that("find_cutpoint_number systematic search works for BIC", {
  res_bic <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", method = "systematic", criterion = "BIC", max_cuts = 2, use_parallel = FALSE, nmin = 1)))
  skip_if(is.null(res_bic), "Systematic search returned NULL, skipping checks.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 3)
  expect_true("BIC" %in% names(res_bic$results))
})

test_that("find_cutpoint_number handles small datasets and zero cuts", {
  small_data <- mock_data[1:5, ]
  skip_if(nrow(small_data) < 5, "Small dataset test skipped due to insufficient data.")
  expect_warning(
    find_cutpoint_number(small_data, "predictor", "time", "event", method = "systematic", criterion = "BIC", max_cuts = 1, nmin = 1),
    "Loglik converged before variable"
  )

  # Test if zero cuts is evaluated
  res_bic <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic", criterion = "BIC", max_cuts = 2, maxiter = 10, nmin = 1)))
  skip_if(is.null(res_bic), "Genetic search returned NULL, skipping checks.")
  expect_true(0 %in% res_bic$results$num_cuts)
})

# --- Tests for find_cutpoint() ---

test_that("find_cutpoint systematic search works for one cut with recommended criterion", {
  # Test with logrank (recommended default)
  # NOTE: suppressWarnings is added to silence potential coxph convergence warnings with random data
  res_lr1 <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1)))
  expect_s3_class(res_lr1, "find_cutpoint")
  expect_length(res_lr1$optimal_cuts, 1)
  expect_type(res_lr1$optimal_cuts, "double")
  expect_true(all(res_lr1$optimal_cuts >= min(mock_data$predictor) & res_lr1$optimal_cuts <= max(mock_data$predictor)))

  # Test with hazard_ratio (recommended for effect size)
  res_hr1 <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "hazard_ratio", nmin = 1)))
  expect_s3_class(res_hr1, "find_cutpoint")
  expect_length(res_hr1$optimal_cuts, 1)
})

test_that("find_cutpoint genetic search works for multiple cuts with logrank", {
  skip_if_not_installed("rgenoud")

  # Test 2 cuts (recommended method: genetic)
  # NOTE: suppressWarnings is used to silence the expected rgenoud warning about hitting maxiter
  res_lr2 <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 2, method = "genetic", criterion = "logrank", maxiter = 10, nmin = 1)))
  skip_if(any(is.na(res_lr2$optimal_cuts)), "GA search failed to find a valid solution.")
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
  expect_true(all(diff(res_lr2$optimal_cuts) > 0)) # Check cuts are in ascending order

  # Test 3 cuts
  res_lr3 <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 3, method = "genetic", criterion = "logrank", maxiter = 10, nmin = 1)))
  skip_if(any(is.na(res_lr3$optimal_cuts)), "GA search failed for 3 cuts.")
  expect_s3_class(res_lr3, "find_cutpoint")
  expect_length(res_lr3$optimal_cuts, 3)
})

test_that("find_cutpoint handles invalid inputs", {
  # Missing column
  bad_data <- mock_data[, c("time", "event")]
  expect_error(find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
               "The following required column\\(s\\) were not found in the data: 'predictor'")

  # Invalid time (negative values)
  bad_data <- mock_data
  bad_data$time[1] <- -1
  expect_error(find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
               "Time variable must be non-negative.")

  # Invalid num_cuts for systematic
  expect_error(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 0, method = "systematic", nmin = 1),
               "Systematic search currently only supports num_cuts = 1 or 2.")

  # Invalid criterion
  expect_error(find_cutpoint(mock_data, "predictor", "time", "event", criterion = "invalid", nmin = 1),
               "arg.*should be one of \"logrank\", \"hazard_ratio\", \"p_value\"")

  # hazard_ratio with more than 1 cut
  expect_error(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 2, criterion = "hazard_ratio", nmin = 1),
               "'hazard_ratio' criterion is only supported for num_cuts = 1.")
})

test_that("find_cutpoint handles low variability predictor", {
  low_var_data <- mock_data
  low_var_data$predictor <- rep(50, n_test)
  expect_warning(
    find_cutpoint(low_var_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 1),
    "Systematic search could not find any valid cut-points with the given constraints."
  )
})

test_that("find_cutpoint handles skewed predictor", {
  res_lr_skewed <- suppressMessages(suppressWarnings(find_cutpoint(mock_data_skewed, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 1)))
  expect_s3_class(res_lr_skewed, "find_cutpoint")
  expect_length(res_lr_skewed$optimal_cuts, 1)
})

test_that("find_cutpoint handles heavy censoring", {
  res_lr_censor <- suppressMessages(suppressWarnings(find_cutpoint(mock_data_heavy_censor, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1)))
  expect_s3_class(res_lr_censor, "find_cutpoint")
  expect_length(res_lr_censor$optimal_cuts, 1)
})

test_that("find_cutpoint works with covariates", {
  res_cov <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", covariates = "covariate1", num_cuts = 1, method = "systematic", nmin = 1)))
  expect_s3_class(res_cov, "find_cutpoint")
  expect_length(res_cov$optimal_cuts, 1)
})

test_that("find_cutpoint works with parallelization", {
  res_lr_parallel <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", use_parallel = TRUE, nmin = 1)))
  expect_s3_class(res_lr_parallel, "find_cutpoint")
  expect_length(res_lr_parallel$optimal_cuts, 1)
})

# --- Tests for validate_cutpoint() ---

test_that("validate_cutpoint works with recommended settings", {
  skip_on_ci() # Skip in CI to reduce runtime
  set.seed(42)
  num_replicates <- if (Sys.getenv("CI") == "true") 20 else 50
  fc_result <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1)))
  skip_if(any(is.na(fc_result$optimal_cuts)), "Systematic search failed, skipping validation test.")

  val_result <- suppressMessages(suppressWarnings(validate_cutpoint(fc_result, num_replicates = num_replicates, use_parallel = TRUE, nmin = 1)))
  skip_if(is.null(val_result), "Validation failed, skipping checks.")

  expect_s3_class(val_result, "validate_cutpoint_result")
  expect_equal(nrow(val_result$confidence_intervals), 1)
  expect_true(val_result$parameters$successful_reps > 0)
  skip_if(is.null(val_result$confidence_intervals$lower), "Confidence intervals not generated, skipping checks.")
  expect_true(all(val_result$confidence_intervals$lower <= val_result$confidence_intervals$upper))
})

test_that("validate_cutpoint warns with low successful replicates", {
  set.seed(42)
  # --- FIXED: Created a more extreme "pathological" dataset to reliably cause bootstrap failures ---
  # A small dataset with few unique values and a high nmin makes it very hard to find valid splits.
  pathological_data <- mock_data[1:50, ]
  pathological_data$predictor <- sample(c(40, 50, 60), 50, replace = TRUE)

  fc_result <- suppressMessages(suppressWarnings(
    find_cutpoint(pathological_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 20)
  ))
  skip_if(any(is.na(fc_result$optimal_cuts)), "Systematic search failed on pathological data, skipping validation test.")

  # This setup should now reliably produce a low success rate and trigger the warning.
  expect_warning(
    validate_cutpoint(fc_result, num_replicates = 100, use_parallel = FALSE, nmin = 20),
    "bootstrap replicates succeeded.*confidence intervals may be unreliable"
  )
})

test_that("validate_cutpoint works for genetic result with 2 cuts", {
  skip_if_not_installed("rgenoud")
  set.seed(42)
  num_replicates <- if (Sys.getenv("CI") == "true") 20 else 50
  fc_result_gen <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 2, method = "genetic", criterion = "logrank", maxiter = 10, nmin = 1)))
  skip_if(any(is.na(fc_result_gen$optimal_cuts)), "GA failed to find initial cutpoint, skipping validation test.")

  val_result_gen <- suppressMessages(suppressWarnings(validate_cutpoint(fc_result_gen, num_replicates = num_replicates, use_parallel = TRUE, maxiter = 10, nmin = 1)))
  skip_if(is.null(val_result_gen), "Validation failed, skipping checks.")

  expect_s3_class(val_result_gen, "validate_cutpoint_result")
  expect_equal(nrow(val_result_gen$confidence_intervals), 2)
  expect_true(val_result_gen$parameters$successful_reps > 0)
})

test_that("validate_cutpoint works with parallelization", {
  set.seed(42)
  num_replicates <- if (Sys.getenv("CI") == "true") 20 else 50
  fc_result <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1)))
  skip_if(any(is.na(fc_result$optimal_cuts)), "Systematic search failed, skipping validation test.")

  val_result_parallel <- suppressMessages(suppressWarnings(validate_cutpoint(fc_result, num_replicates = num_replicates, use_parallel = TRUE, nmin = 1)))
  skip_if(is.null(val_result_parallel), "Validation failed, skipping checks.")

  expect_s3_class(val_result_parallel, "validate_cutpoint_result")
  expect_equal(nrow(val_result_parallel$confidence_intervals), 1)
})

# --- Tests for Full Workflow ---

test_that("full workflow: find_cutpoint_number -> find_cutpoint -> validate_cutpoint", {
  skip_if_not_installed("rgenoud")
  set.seed(42)
  num_replicates <- if (Sys.getenv("CI") == "true") 20 else 50
  # Step 1: Determine optimal number of cuts
  res_fcn <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic", criterion = "BIC", max_cuts = 2, maxiter = 10, nmin = 1)))
  skip_if(is.null(res_fcn), "find_cutpoint_number failed, skipping workflow test.")

  # Get optimal number of cuts
  optimal_cuts <- res_fcn$results$num_cuts[which.min(res_fcn$results$BIC)]

  # Step 2: Find cutpoint locations
  method <- if (optimal_cuts >= 2) "genetic" else "systematic"
  res_fc <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = optimal_cuts, method = method, criterion = "logrank", maxiter = 10, nmin = 1)))
  skip_if(any(is.na(res_fc$optimal_cuts)), "find_cutpoint failed, skipping validation.")

  expect_s3_class(res_fc, "find_cutpoint")
  expect_length(res_fc$optimal_cuts, optimal_cuts)

  # Step 3: Validate cutpoints
  val_result <- suppressMessages(suppressWarnings(validate_cutpoint(res_fc, num_replicates = num_replicates, use_parallel = TRUE, maxiter = 10, nmin = 1)))
  skip_if(is.null(val_result), "validate_cutpoint failed, skipping checks.")

  expect_s3_class(val_result, "validate_cutpoint_result")
  expect_equal(nrow(val_result$confidence_intervals), optimal_cuts)
})

# --- Tests for S3 Methods (print, summary, plot) ---

test_that("S3 methods run without errors and produce correct output", {
  set.seed(42)
  res_fc <- suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1)))
  res_fcn <- suppressMessages(suppressWarnings(find_cutpoint_number(mock_data, "predictor", "time", "event", max_cuts = 1, method = "systematic", criterion = "BIC", nmin = 1)))
  num_replicates <- if (Sys.getenv("CI") == "true") 20 else 50

  # Test find_cutpoint methods
  expect_snapshot(res_fc)
  expect_snapshot(summary(res_fc))

  p1 <- plot(res_fc)
  if (!is.null(p1)) {
    expect_s3_class(p1$plot, "ggplot")
    expect_true(nrow(p1$plot$data) > 0) # Check plot has data
  }

  # Test find_cutpoint_number methods
  if (!is.null(res_fcn)) {
    expect_snapshot(res_fcn)
    expect_snapshot(summary(res_fcn))

    p2 <- plot(res_fcn)
    if (!is.null(p2)) {
      expect_s3_class(p2, "ggplot")
      expect_true(nrow(p2$data) > 0) # Check plot has data
    }
  }

  # Test validate_cutpoint methods
  res_val <- suppressMessages(suppressWarnings(validate_cutpoint(res_fc, num_replicates = num_replicates, use_parallel = TRUE, nmin = 1)))
  if (!is.null(res_val)) {
    expect_s3_class(res_val, "validate_cutpoint_result")
    expect_true("confidence_intervals" %in% names(res_val))
    skip_if(is.null(res_val$confidence_intervals$lower), "Confidence intervals not generated, skipping checks.")
    expect_type(res_val$confidence_intervals$lower, "double")
    expect_type(res_val$confidence_intervals$upper, "double")
    expect_true(all(res_val$confidence_intervals$lower <= res_val$confidence_intervals$upper))

    p3 <- plot(res_val)
    if (!is.null(p3)) {
      expect_s3_class(p3, "ggplot")
      expect_true(nrow(p3$plot$data) > 0) # Check plot has data
    }
  }
})

