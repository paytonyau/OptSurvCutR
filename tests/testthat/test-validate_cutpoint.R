# ===================================================================
# TESTS: validate_cutpoint()
# Bootstrap resampling, stability assessment, and parallel processing
# ===================================================================

# --- Core Validation & Bootstrapping ---

#' @srrstats {RE7.0} Tests bootstrap resampling.
#' @srrstats {RE7.1} Tests cut-point stability assessment.
test_that("validate_cutpoint works with recommended settings", {
  # Capture messages to verify the stability message appears
  expect_message(
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot,
                                  num_replicates = 20,
                                  n_cores = 1,
                                  max.generations = 10,
                                  seed = 42
      )
    }),
    regexp = "90% of original"
  )

  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(result$parameters$successful_reps >= 10)
})

#' @srrstats {RE7.0} Tests bootstrap resampling.
#' @srrstats {G5.8} Tests validation of a stochastic algorithm result.
test_that("validate_cutpoint works for genetic result with 2 cuts", {
  skip_if_not_installed("rgenoud")
  fc_result_gen2 <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data_3groups,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 2,
      method = "genetic",
      nmin = 5,
      max.generations = 5,
      quiet = TRUE
    )
  ))
  skip_if(
    any(is.na(fc_result_gen2$optimal_cuts)),
    "GA 2-cut search failed, skipping validation test."
  )
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(fc_result_gen2,
                                  num_replicates = 20,
                                  n_cores = 1,
                                  nmin = 5,
                                  max.generations = 10,
                                  seed = 42
      )
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(ncol(result$bootstrap_distribution) == 2)
  expect_true(result$parameters$successful_reps >= 10)
})

#' @srrstats {G5.6} Tests optional parallel via doParallel.
test_that("validate_cutpoint works with parallelisation", {
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot,
                                  num_replicates = 20,
                                  n_cores = 2,
                                  nmin = 2,
                                  max.generations = 10,
                                  seed = 42
      )
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(result$parameters$successful_reps >= 10)
})

#' @srrstats {RE7.0} Tests bootstrap resampling.
#' @srrstats {RE2.2} Tests covariate adjustment within bootstrap.
test_that("validate_cutpoint works with covariates", {
  fc_result_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      covariates = "covariate1",
      nmin = 5,
      quiet = TRUE
    )
  ))
  skip_if(
    any(is.na(fc_result_cov$optimal_cuts)),
    "find_cutpoint with covariates failed."
  )
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(fc_result_cov,
                                  num_replicates = 20,
                                  n_cores = 1,
                                  nmin = 5,
                                  max.generations = 10,
                                  seed = 42
      )
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(result$parameters$successful_reps >= 10)
})

# --- Input Validation & Edge Cases ---

#' @srrstats {G2.0} Tests input object class validation.
#' @srrstats {G2.1a} Tests check for valid cutpoint_result structure.
#' @srrstats {G2.2} Tests that NA in optimal_cuts triggers error.
#' @srrstats {G2.13} Tests cli_abort for invalid input.
test_that("validate_cutpoint handles invalid inputs", {
  expect_error(
    validate_cutpoint(list(), num_replicates = 50, nmin = 1),
    "Input must be a `find_cutpoint` object"
  )

  fc_result_na <- valid_fc_result_for_boot
  fc_result_na$optimal_cuts <- NA
  expect_error(
    validate_cutpoint(fc_result_na, num_replicates = 50, nmin = 1),
    regexp = "Input `find_cutpoint` object has NA cut-points"
  )
})

#' @srrstats {G2.0a} Tests validation of scalar parameter num_replicates.
#' @srrstats {G5.2a} Tests handling of <20 successful replicates.
test_that("validate_cutpoint handles low replicates", {
  expect_error(
    validate_cutpoint(valid_fc_result_for_boot,
                      num_replicates = 5,
                      n_cores = 1, nmin = 5
    ),
    regexp = "`num_replicates` must be >= 20 for validation\\."
  )
})

#' @srrstats {G2.0a} Tests validation of scalar parameter num_replicates.
test_that("Coverage: validate_cutpoint - non-integer num_replicates", {
  expect_error(
    validate_cutpoint(valid_fc_result_for_boot, num_replicates = 25.5, n_cores = 1),
    regexp = "positive integer"
  )
})

#' @srrstats {G2.0a} Tests validation of scalar parameter num_replicates.
#' @srrstats {G5.8d} Tests data outside scope (insufficient sample size for requested cuts).
test_that("Coverage: validate_cutpoint - insufficient data after nmin", {
  valid_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data[1:30, ], "predictor", "time", "event",
                  num_cuts = 1, nmin = 5, quiet = TRUE
    )
  ))
  skip_if(any(is.na(valid_fc$optimal_cuts)))
  expect_error(
    validate_cutpoint(valid_fc,
                      num_replicates = 20,
                      n_cores = 1,
                      nmin = 16
    ),
    regexp = "Not enough data"
  )
})

test_that("Coverage: validate_cutpoint - parallel with n_cores < 1", {
  result <- suppressMessages(suppressWarnings(
    validate_cutpoint(valid_fc_result_for_boot,
                      num_replicates = 20,
                      n_cores = 0, nmin = 5, seed = 42
    )
  ))
  expect_equal(result$parameters$n_cores, 1)
})

test_that("Coverage: validate_cutpoint failure modes (Pathological Data)", {
  fc_result_pathological <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data_pathological,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      nmin = 5,
      quiet = TRUE
    )
  ))
  expect_true(is.na(fc_result_pathological$optimal_cuts[1]))

  expect_error(
    validate_cutpoint(
      fc_result_pathological,
      num_replicates = 20, seed = 1, n_cores = 1
    ),
    regexp = "Input `find_cutpoint` object has NA cut-points"
  )
})

# --- S3 Methods for validate_cutpoint ---

#' @srrstats {RE4.17} Tests S3 methods (print, summary, plot).
#' @srrstats {RE6.0} Tests diagnostic plot of bootstrap distribution.
test_that("S3 methods for validate_cutpoint handle arguments", {
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot,
                                  num_replicates = 20,
                                  n_cores = 1,
                                  nmin = 5,
                                  max.generations = 10,
                                  seed = 42
      )
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")

  expect_output(
    print(result),
    regexp = paste(
      "Cut-point Stability Analysis",
      "Successful Replicates",
      "Confidence Intervals",
      "Bootstrap Summary Statistics",
      sep = ".*"
    )
  )

  expect_output(
    summary(
      result,
      show_descriptives = TRUE, show_ci = TRUE, show_params = TRUE
    ),
    regexp = paste(
      "Cut-point Stability Analysis",
      "Bootstrap Distribution Summary",
      "Confidence Intervals",
      "Validation Parameters",
      sep = ".*"
    )
  )

  expect_s3_class(plot(result), "ggplot")
})

test_that("Coverage: validate_cutpoint S3 methods for 0 successes", {
  res_val_empty <- structure(
    list(
      bootstrap_distribution = data.frame(),
      parameters = list(successful_reps = 0, num_replicates = 20),
      original_cuts = 50,
      confidence_intervals = data.frame(Lower = NA, Upper = NA),
      boot_summary = list(Cut1 = list(
        mean = NA, sd = NA, median = NA, Q1 = NA, Q3 = NA
      ))
    ),
    class = "validate_cutpoint_result"
  )

  expect_message(
    plot(res_val_empty),
    regexp = "Cannot plot: 0 successful bootstrap replicates"
  )

  expect_output(
    summary(res_val_empty, show_descriptives = FALSE),
    "Successful Replicates: 0 / 20"
  )
})
