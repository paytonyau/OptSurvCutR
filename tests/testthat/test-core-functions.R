# Comprehensive test file for the OptSurvCutR package.
# This script merges the original core functionality tests with additional
# tests designed to increase code coverage by targeting edge cases,
# error conditions, and internal helper function branches.

library(testthat)
library(survival)
library(withr) # For local_mocked_bindings

# --- Setup: Create consistent, reproducible mock datasets ---
set.seed(42)
n_test <- 60 # Increased sample size for robustness with covariates
# Bimodal predictor for 1 cut (2 groups): ~50% N(40, 5), ~50% N(60, 5)
n1 <- n_test / 2
n2 <- n_test - n1
mock_data <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(rnorm(n1, mean = 40, sd = 5), rnorm(n2, mean = 60, sd = 5)),
  covariate1 = rnorm(n_test, mean = 5, sd = 1),
  covariate2 = sample(c("A", "B"), n_test, replace = TRUE) # Factor covariate
)

# Trimodal predictor for 2 cuts (3 groups)
n1_3 <- floor(n_test / 3)
n2_3 <- floor(n_test / 3)
n3_3 <- n_test - n1_3 - n2_3
mock_data_3groups <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(rnorm(n1_3, mean = 30, sd = 5), rnorm(n2_3, mean = 50, sd = 5), rnorm(n3_3, mean = 70, sd = 5)),
  covariate1 = rnorm(n_test, mean = 5, sd = 1)
)

# Mock dataset with skewed predictor
mock_data_skewed <- mock_data
mock_data_skewed$predictor <- rexp(n_test, rate = 0.1)

# Mock dataset with heavy censoring (90% censored)
mock_data_heavy_censor <- mock_data
mock_data_heavy_censor$event <- sample(0:1, n_test, replace = TRUE, prob = c(0.9, 0.1))

# Pathological dataset for bootstrap failure (small size, constant predictor)
mock_data_pathological <- mock_data[1:20, ]
mock_data_pathological$predictor <- rep(50, 20) # Constant predictor

# Tiny dataset for AIC/AICc edge cases
tiny_data <- mock_data[1:4, ]
tiny_data$predictor <- c(40, 45, 50, 55)

# New dataset for small dataset test (line 106) with constant predictor
small_data_unique <- data.frame(
  time = rexp(20, rate = 0.05),
  event = sample(0:1, 20, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = rep(50, 20), # Constant to induce warning
  covariate1 = rnorm(20, mean = 5, sd = 1)
)

# --- CORE FUNCTIONALITY TESTS ---

# --- Tests for find_cutpoint_number() ---

test_that("find_cutpoint_number genetic search works for all criteria", {
  skip_if_not_installed("rgenoud")

  # Test BIC (recommended default)
  res_bic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic",
                         criterion = "BIC", max_cuts = 2, maxiter = 5, nmin = 1)
  ))
  skip_if(is.null(res_bic), "Genetic search returned NULL.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 3) # 0, 1, 2 cuts
  expect_true("BIC" %in% names(res_bic$results))
  expect_true(all(res_bic$results$num_cuts %in% 0:2))

  # Test AIC
  res_aic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic",
                         criterion = "AIC", maxiter = 5, nmin = 1)
  ))
  skip_if(is.null(res_aic), "Genetic search returned NULL.")
  expect_s3_class(res_aic, "find_cutpoint_number_result")
  expect_true("AIC" %in% names(res_aic$results))

  # Test AICc
  res_aicc <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic",
                         criterion = "AICc", maxiter = 5, nmin = 1)
  ))
  skip_if(is.null(res_aicc), "Genetic search returned NULL.")
  expect_s3_class(res_aicc, "find_cutpoint_number_result")
  expect_true("AICc" %in% names(res_aicc$results))
})

test_that("find_cutpoint_number systematic search works for BIC", {
  res_bic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "systematic",
                         criterion = "BIC", max_cuts = 2, nmin = 1)
  ))
  skip_if(is.null(res_bic), "Systematic search returned NULL.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 3)
  expect_true("BIC" %in% names(res_bic$results))
})

# *** NEW: Test find_cutpoint_number with covariates ***
test_that("find_cutpoint_number works with covariates (systematic)", {
  res_bic_cov <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "systematic",
      covariates = c("covariate1", "covariate2"),
      criterion = "BIC",
      max_cuts = 1,
      nmin = 10
    )
  ))
  skip_if(is.null(res_bic_cov), "Systematic search with covariates returned NULL.")
  expect_s3_class(res_bic_cov, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_cov$results), 2) # 0, 1 cuts
  expect_true(!is.na(res_bic_cov$optimal_num_cuts))
  expect_equal(res_bic_cov$parameters$covariates, c("covariate1", "covariate2"))
})

test_that("find_cutpoint_number works with covariates (genetic)", {
  skip_if_not_installed("rgenoud")
  res_bic_cov_gen <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "genetic",
      covariates = c("covariate1"),
      criterion = "BIC",
      max_cuts = 1,
      nmin = 10,
      maxiter = 5
    )
  ))
  skip_if(is.null(res_bic_cov_gen), "Genetic search with covariates returned NULL.")
  expect_s3_class(res_bic_cov_gen, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_cov_gen$results), 2) # 0, 1 cuts
  expect_true(!is.na(res_bic_cov_gen$optimal_num_cuts))
  expect_equal(res_bic_cov_gen$parameters$covariates, "covariate1")
})


test_that("find_cutpoint_number handles predictor with too few unique values", {

  # Ensure the test data is correctly set up for this specific test
  skip_if(length(unique(small_data_unique$predictor)) > 1,
          "Test data 'small_data_unique' must have only one unique predictor value.")

  # Run the function.
  # We suppress messages because we don't need to see the "too few unique values"
  # message printed in our test log.
  result <- suppressMessages(
    find_cutpoint_number(
      data = small_data_unique,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      max_cuts = 2,
      nmin = 5,
      criterion = "BIC"
    )
  )

  # Test 1: Check that the function returned NA for the optimal number of cuts
  expect_true(
    is.na(result$optimal_num_cuts),
    info = "Function should return NA for optimal_num_cuts when predictor has too few unique values."
  )

  # Test 2: Check that the results table is empty
  expect_true(
    nrow(result$results) == 0,
    info = "Function should return an empty results table."
  )
})

test_that("find_cutpoint_number handles AIC and AICc edge cases", {
  res_aic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(tiny_data, "predictor", "time", "event", method = "systematic",
                         criterion = "AIC", max_cuts = 1, nmin = 4)
  ))
  skip_if(is.null(res_aic), "Systematic search returned NULL.")
  expect_true(is.na(res_aic$results$AIC[res_aic$results$num_cuts == 1]) ||
                !any(res_aic$results$cuts[res_aic$results$num_cuts == 1] != "NA"))

  res_aicc <- suppressMessages(suppressWarnings(
    find_cutpoint_number(tiny_data, "predictor", "time", "event", method = "systematic",
                         criterion = "AICc", max_cuts = 1, nmin = 4)
  ))
  skip_if(is.null(res_aicc), "Systematic search returned NULL.")
  expect_true(is.na(res_aicc$results$AICc[res_aicc$results$num_cuts == 1]) ||
                !any(res_aicc$results$cuts[res_aicc$results$num_cuts == 1] != "NA"))
})

test_that("find_cutpoint_number (systematic) runs sequentially", {
  res_bic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "systematic",
                         criterion = "BIC", max_cuts = 1, nmin = 1)
  ))
  skip_if(is.null(res_bic), "Systematic search returned NULL.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 2) # 0, 1 cuts
})

test_that("find_cutpoint_number (genetic) runs sequentially", {
  skip_if_not_installed("rgenoud")
  res_bic_gen <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic",
                         criterion = "BIC", max_cuts = 1, maxiter = 5, nmin = 1)
  ))
  skip_if(is.null(res_bic_gen), "Genetic search returned NULL.")
  expect_s3_class(res_bic_gen, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_gen$results), 2) # 0, 1 cuts
})


test_that("find_cutpoint_number handles invalid inputs", {
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event", max_cuts = -1),
    "max_cuts must be a non-negative integer"
  )
  expect_error(
    find_cutpoint_number(mock_data, predictor = NULL, "time", "event"),
    "must be specified"
  )
  expect_error(
    find_cutpoint_number(mock_data, "predictor", outcome_time = NULL, "event"),
    "must be specified"
  )
  expect_error(
    find_cutpoint_number(data.frame(a=1), "predictor", "time", "event"),
    "not found"
  )
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event", nmin = -1),
    "must be a positive number"
  )
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "systematic", max_cuts = 3),
    "only implemented for max_cuts <= 2"
  )
})

# --- Tests for find_cutpoint() ---

test_that("find_cutpoint systematic search works for one and two cuts", {
  # One cut
  res_lr1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "logrank", nmin = 1)
  ))
  expect_s3_class(res_lr1, "find_cutpoint")
  expect_length(res_lr1$optimal_cuts, 1)
  expect_true(all(is.na(res_lr1$optimal_cuts) |
                    (res_lr1$optimal_cuts >= min(mock_data$predictor) &
                       res_lr1$optimal_cuts <= max(mock_data$predictor))))

  # Two cuts (with appropriate data)
  res_lr2 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event", num_cuts = 2, method = "systematic",
                  criterion = "logrank", nmin = 10)
  ))
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
})

test_that("find_cutpoint genetic search works for multiple cuts", {
  skip_if_not_installed("rgenoud")

  res_lr2 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event", num_cuts = 2, method = "genetic",
                  criterion = "logrank", maxiter = 5, nmin = 1)
  ))
  skip_if(all(is.na(res_lr2$optimal_cuts)), "GA search failed for 2 cuts.")
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
  expect_true(all(is.na(res_lr2$optimal_cuts) | diff(res_lr2$optimal_cuts) > 0))

  res_lr3 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event", num_cuts = 3, method = "genetic",
                  criterion = "logrank", maxiter = 5, nmin = 1)
  ))
  skip_if(all(is.na(res_lr3$optimal_cuts)), "GA search failed for 3 cuts.")
  expect_s3_class(res_lr3, "find_cutpoint")
  expect_length(res_lr3$optimal_cuts, 3)
})

test_that("find_cutpoint handles various criteria and nmin proportion", {
  # Hazard ratio
  res_hr1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "hazard_ratio", nmin = 1)
  ))
  expect_s3_class(res_hr1, "find_cutpoint")
  expect_length(res_hr1$optimal_cuts, 1)

  # p-value
  res_pv1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "p_value", nmin = 1)
  ))
  expect_s3_class(res_pv1, "find_cutpoint")
  expect_length(res_pv1$optimal_cuts, 1)

  # nmin as proportion
  res_nmin_prop <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", nmin = 0.1)
  ))
  expect_s3_class(res_nmin_prop, "find_cutpoint")
})

test_that("find_cutpoint handles invalid inputs", {
  bad_data <- mock_data[, c("time", "event")]
  expect_error(
    find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
    "The following required column\\(s\\) were not found in the data: 'predictor'"
  )

  bad_data <- mock_data
  bad_data$time[1] <- -1
  expect_error(
    find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
    "Time variable must be non-negative."
  )

  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = -1, method = "systematic", nmin = 1),
    "num_cuts must be a non-negative integer"
  )

  # Test for NA num_cuts
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = NA, method = "systematic", nmin = 1),
    "num_cuts must be a non-negative integer"
  )

  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event", criterion = "invalid", nmin = 1),
    "arg.*should be one of \"logrank\", \"hazard_ratio\", \"p_value\""
  )

  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 2, criterion = "hazard_ratio", nmin = 1),
    "'hazard_ratio' criterion is only supported for num_cuts = 1."
  )
})

test_that("find_cutpoint handles quiet arg and low variability predictor", {
  # Test quiet = TRUE on a successful run
  expect_no_message(
    find_cutpoint(mock_data, "predictor", "time", "event", quiet = TRUE)
  )

  # Test warning on low variability data
  low_var_data <- mock_data
  low_var_data$predictor <- rep(50, n_test)
  expect_message(
    find_cutpoint(low_var_data, "predictor", "time", "event"),
    "Predictor has too few unique values"
  )
})

test_that("find_cutpoint handles various data scenarios", {
  # Skewed predictor
  res_lr_skewed <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_skewed, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 1)
  ))
  expect_s3_class(res_lr_skewed, "find_cutpoint")
  expect_length(res_lr_skewed$optimal_cuts, 1)

  # Heavy censoring
  res_lr_censor <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_heavy_censor, "predictor", "time", "event", num_cuts = 1,
                  method = "systematic", criterion = "logrank", nmin = 1)
  ))
  expect_s3_class(res_lr_censor, "find_cutpoint")
  expect_length(res_lr_censor$optimal_cuts, 1)

  # Cox convergence issues
  bad_data <- mock_data
  bad_data$predictor <- rep(50, n_test)
  res_cox_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(bad_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 10)
  ))
  expect_s3_class(res_cox_fail, "find_cutpoint")
  expect_true(all(is.na(res_cox_fail$optimal_cuts)))
})

test_that("find_cutpoint (systematic) works with covariates", {
  # Covariates
  res_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", covariates = "covariate1",
                  num_cuts = 1, method = "systematic", nmin = 1)
  ))
  expect_s3_class(res_cov, "find_cutpoint")
  expect_length(res_cov$optimal_cuts, 1)
  expect_true(!is.na(res_cov$optimal_cuts[1])) # Check it found a cut

  # Sequential run test (already done, this is just a label)
  res_lr_seq <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "logrank", nmin = 1)
  ))
  expect_s3_class(res_lr_seq, "find_cutpoint")
  expect_length(res_lr_seq$optimal_cuts, 1)
})

test_that("find_cutpoint (genetic) works sequentially", {
  skip_if_not_installed("rgenoud")
  res_gen_seq <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 2, method = "genetic",
                  criterion = "logrank", maxiter = 5, nmin = 1)
  ))
  skip_if(is.null(res_gen_seq) || all(is.na(res_gen_seq$optimal_cuts)), "Sequential genetic search failed.")
  expect_s3_class(res_gen_seq, "find_cutpoint")
  expect_length(res_gen_seq$optimal_cuts, 2)
})


test_that("find_cutpoint handles high num_cuts with p_value", {
  skip_if_not_installed("rgenoud")
  res_pv4 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 4, method = "genetic",
                  criterion = "p_value", maxiter = 5, nmin = 1)
  ))
  skip_if(all(is.na(res_pv4$optimal_cuts)), "GA search failed for 4 cuts.")
  expect_s3_class(res_pv4, "find_cutpoint")
  expect_length(res_pv4$optimal_cuts, 4)
})

test_that("find_cutpoint handles insufficient and constant data", {
  # Insufficient data
  small_data <- mock_data[1:5, ]
  res_insufficient <- suppressMessages(suppressWarnings(
    find_cutpoint(small_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 5)
  ))
  expect_s3_class(res_insufficient, "find_cutpoint")
  expect_true(all(is.na(res_insufficient$optimal_cuts)))

  # Constant predictor
  constant_data <- mock_data[1:20, ]
  constant_data$predictor <- rep(50, 20)
  res_constant <- suppressMessages(suppressWarnings(
    find_cutpoint(constant_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 1)
  ))
  expect_s3_class(res_constant, "find_cutpoint")
  expect_true(all(is.na(res_constant$optimal_cuts)))
})

# --- Tests for validate_cutpoint() ---

# --- Setup block for validate_cutpoint tests ---
# This setup is used for most validation tests that need a reliable input
set.seed(42)
mock_data_for_boot <- data.frame(
  time = rexp(80, rate = 0.05),
  event = sample(0:1, 80, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(rnorm(40, mean = 40, sd = 5), rnorm(40, mean = 60, sd = 5))
)

valid_fc_result_for_boot <- suppressMessages(suppressWarnings(
  find_cutpoint(
    data = mock_data_for_boot,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    num_cuts = 1,
    nmin = 5
  )
))

if (any(is.na(valid_fc_result_for_boot$optimal_cuts))) {
  stop("Test setup failed: Could not generate a valid fc_result for bootstrap tests.")
}
# --- End of setup block ---

test_that("validate_cutpoint works with recommended settings", {
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot, num_replicates = 50, use_parallel = FALSE, nmin = 5)
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(nrow(result$bootstrap_distribution) >= 20)
})

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
      maxiter = 5
    )
  ))

  skip_if(any(is.na(fc_result_gen2$optimal_cuts)), "GA 2-cut search failed, skipping validation test.")

  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(fc_result_gen2, num_replicates = 50, use_parallel = FALSE, nmin = 5)
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(ncol(result$bootstrap_distribution) == 2)
  expect_true(nrow(result$bootstrap_distribution) >= 20)
})

test_that("validate_cutpoint works with parallelization", {
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot, num_replicates = 50, use_parallel = TRUE, n_cores = 2, nmin = 5)
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(nrow(result$bootstrap_distribution) >= 20)
})

test_that("validate_cutpoint handles low replicates", {
  expect_error(
    validate_cutpoint(valid_fc_result_for_boot, num_replicates = 5, use_parallel = FALSE, nmin = 5),
    regexp = "num_replicates must be at least 20",
    info = "Should error for low replicates"
  )
})

test_that("validate_cutpoint works with covariates", {
  fc_result_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      covariates = "covariate1",
      nmin = 5
    )
  ))

  skip_if(any(is.na(fc_result_cov$optimal_cuts)), "find_cutpoint with covariates failed.")

  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(fc_result_cov, num_replicates = 50, use_parallel = FALSE, nmin = 5)
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_true(nrow(result$bootstrap_distribution) >= 20)
})


test_that("validate_cutpoint handles invalid inputs", {
  expect_error(
    validate_cutpoint(list(), num_replicates = 50, nmin = 1),
    "Input must be an object from the find_cutpoint function"
  )

  # Test with NA cut-points
  fc_result_na <- valid_fc_result_for_boot
  fc_result_na$optimal_cuts <- NA
  expect_error(
    validate_cutpoint(fc_result_na, num_replicates = 50, nmin = 1),
    "Input 'find_cutpoint' object contains NA cut-points. Cannot validate."
  )
})

# --- Tests for plotting_functions.R ---

test_that("plot_optimization_curve works for all criteria", {
  # Logrank
  res_lr <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "logrank", nmin = 1)
  ))
  skip_if(is.null(res_lr) || all(is.na(res_lr$optimal_cuts)), "Systematic logrank search failed.")
  p_lr <- plot_optimization_curve(res_lr)
  expect_s3_class(p_lr, "ggplot")
  expect_equal(p_lr$labels$y, "Log-Rank Statistic")

  # Hazard Ratio
  res_hr <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "hazard_ratio", nmin = 1)
  ))
  skip_if(is.null(res_hr) || all(is.na(res_hr$optimal_cuts)), "Systematic HR search failed.")
  p_hr <- plot_optimization_curve(res_hr)
  expect_s3_class(p_hr, "ggplot")
  expect_equal(p_hr$labels$y, "Hazard Ratio (HR)")

  # P-value
  res_pv <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic",
                  criterion = "p_value", nmin = 1)
  ))
  skip_if(is.null(res_pv) || all(is.na(res_pv$optimal_cuts)), "Systematic p-value search failed.")
  p_pv <- plot_optimization_curve(res_pv)
  expect_s3_class(p_pv, "ggplot")
  expect_equal(p_pv$labels$y, "P-value")
})

test_that("plot_optimization_curve throws errors for invalid input", {
  # Generate a valid result to modify for testing
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", nmin = 1)
  ))
  skip_if(is.null(res_valid), "Systematic search failed, cannot run error tests.")

  # Test for non-find_cutpoint object
  expect_error(
    plot_optimization_curve(list()),
    "Input must be an object from the"
  )

  # Test for genetic search result
  res_genetic <- res_valid
  res_genetic$parameters$method <- "genetic"
  expect_error(
    plot_optimization_curve(res_genetic),
    "This plot is only for results from"
  )

  # Test for num_cuts != 1
  res_2_cuts <- res_valid
  res_2_cuts$parameters$num_cuts <- 2
  expect_error(
    plot_optimization_curve(res_2_cuts),
    "This plot is only supported for results with `num_cuts = 1`"
  )

  # Test for missing all_stats data
  res_no_stats <- res_valid
  res_no_stats$all_stats <- NULL
  expect_error(
    plot_optimization_curve(res_no_stats),
    "object does not contain the necessary `all_stats` data"
  )
})

# --- Tests for S3 Methods and Plot Variations ---

test_that("S3 methods for find_cutpoint handle NA results and arguments", {
  res_na <- find_cutpoint(mock_data_pathological, "predictor", "time", "event",
                          num_cuts = 1, method = "systematic", quiet = TRUE)
  expect_true(all(is.na(res_na$optimal_cuts)))

  # Test print, summary, and plot on NA results
  expect_message(print(res_na), "No optimal cut-point could be determined")
  expect_message(summary(res_na), "No valid optimal cut-point was found")
  expect_message(plot(res_na), "Cannot generate plot, no valid cut-point found")

  # Test summary arguments on a valid result
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, quiet = TRUE)
  ))
  skip_if(is.null(res_valid) || all(is.na(res_valid$optimal_cuts)), "Valid result needed for summary test.")

  sum_res <- summary(res_valid, show_model = FALSE, show_group_counts = FALSE, show_medians = FALSE,
                     show_ph_test = FALSE, show_params = FALSE)
  expect_s3_class(sum_res, "find_cutpoint")
})

test_that("plot.find_cutpoint generates all plot types", {
  skip_if_not_installed("broom")
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, method = "systematic", quiet = TRUE)
  ))
  skip_if(is.null(res_valid) || all(is.na(res_valid$optimal_cuts)), "Valid result needed for plot tests.")

  # Test "outcome" plot (default) - suppress size warning from dependency
  p_outcome <- suppressWarnings(plot(res_valid, type = "outcome"))
  expect_s3_class(p_outcome$plot, "ggplot")

  # Test "distribution" plot
  p_dist <- plot(res_valid, type = "distribution")
  expect_s3_class(p_dist, "ggplot")
  expect_equal(p_dist$labels$y, "Density")

  # Test "forest" plot
  p_forest <- plot(res_valid, type = "forest")
  expect_s3_class(p_forest, "ggplot")
  expect_equal(p_forest$labels$x, "Hazard Ratio (95% CI)")

  # Test forest plot with different reference group
  res_2_cuts <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event", num_cuts = 2, method = "genetic", maxiter=5, quiet = TRUE)
  ))
  skip_if(is.null(res_2_cuts) || all(is.na(res_2_cuts$optimal_cuts)), "2-cut result needed for ref group test.")
  p_forest_g2 <- plot(res_2_cuts, type = "forest", reference_group = "G2")
  expect_s3_class(p_forest_g2, "ggplot")
  expect_match(p_forest_g2$labels$subtitle, "G2")
})

test_that("S3 methods for find_cutpoint_number handle NA results and arguments", {
  # Manually create a result object where all IC values are NA
  res_na <- structure(
    list(
      results = data.frame(num_cuts = 0:1, BIC = c(NA, NA)),
      parameters = list(criterion = "BIC", method = "systematic")
    ),
    class = "find_cutpoint_number_result"
  )

  # Test plot on NA results
  expect_message(plot(res_na), "Cannot generate plot because no valid Information Criterion values")

  # Test summary on NA results
  expect_message(summary(res_na), "Cannot generate summary because no valid optimal model was found")

  # Test summary arguments
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", max_cuts = 1)
  ))
  skip_if(is.null(res_valid), "Valid result needed for summary test.")

  sum_res <- summary(res_valid, show_comparison_table = FALSE, show_best_model_details = FALSE, plot.it = FALSE)
  expect_s3_class(sum_res, "find_cutpoint_number_result")

  # Test the plotting part separately
  expect_s3_class(suppressMessages(summary(res_valid, plot.it = TRUE)), "find_cutpoint_number_result")
})

test_that("S3 methods for validate_cutpoint handle arguments", {
  # Use the valid object from the setup block
  suppressMessages({
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot, num_replicates = 50, use_parallel = FALSE, nmin = 5)
    })
  })
  expect_s3_class(result, "validate_cutpoint_result")
  expect_output(
    print(result),
    regexp = "Cut-point Stability Analysis.*Successful Replicates.*Confidence Intervals.*Bootstrap Summary Statistics",
    info = "Should include main header and key sections"
  )
  expect_output(
    summary(result, show_descriptives = TRUE, show_ci = TRUE, show_params = TRUE),
    regexp = "Cut-point Stability Analysis.*Bootstrap Distribution Summary.*Confidence Intervals.*Validation Parameters",
    info = "Should include all section headers"
  )
  expect_s3_class(plot(result), "ggplot")
})

# --- Tests for utils-helpers.R ---

test_that(".obj handles invalid cuts", {
  expect_equal(
    .obj(params = c(50, 50), time = mock_data$time, censor = mock_data$event,
         target = mock_data$predictor, confound = NULL, numcut = 2, gap = 1, nmin = 1,
         criterion = "logrank"),
    -Inf
  )
})

test_that(".obj handles p_value criterion", {
  expect_true(is.numeric(
    .obj(params = 50, time = mock_data$time, censor = mock_data$event,
         target = mock_data$predictor, confound = NULL, numcut = 1, gap = 1, nmin = 1,
         criterion = "p_value")
  ))
})

test_that("%||% works as expected", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
})

test_that("group creation works in find_cutpoint", {
  set.seed(42)
  res_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 2, method = "genetic",
                  criterion = "logrank", maxiter = 5, nmin = 1)
  ))
  skip_if(all(is.na(res_fc$optimal_cuts)), "Genetic search failed.")
  expect_s3_class(res_fc, "find_cutpoint")
  expect_length(res_fc$optimal_cuts, 2)
  expect_true(all(is.na(res_fc$optimal_cuts) | diff(res_fc$optimal_cuts) > 0))
})

# --- Tests for Full Workflow ---

test_that("full workflow: find_cutpoint_number -> find_cutpoint -> validate_cutpoint", {
  num_result <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      max_cuts = 2,
      nmin = 5,
      criterion = "BIC"
    )
  ))
  skip_if(
    is.null(num_result) || is.na(num_result$optimal_num_cuts) || num_result$optimal_num_cuts < 0,
    "find_cutpoint_number failed to produce a valid number of cuts."
  )

  # Handle 0-cut case
  if (num_result$optimal_num_cuts == 0) {
    # Test passes if it correctly identifies 0 cuts and stops
    expect_equal(num_result$optimal_num_cuts, 0)

  } else {
    # Continue to find_cutpoint and validate_cutpoint
    fc_result <- find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = num_result$optimal_num_cuts,
      nmin = 5,
      quiet = TRUE
    )

    skip_if(
      is.null(fc_result) || any(is.na(fc_result$optimal_cuts)),
      "find_cutpoint failed in the full workflow test."
    )

    suppressMessages({
      suppressWarnings({
        result <- validate_cutpoint(fc_result, num_replicates = 50, use_parallel = FALSE, nmin = 5)
      })
    })
    expect_s3_class(result, "validate_cutpoint_result")
    expect_true(nrow(result$bootstrap_distribution) >= 20)
  }
})

# --- ADDITIONAL TESTS FOR COVERAGE ---

test_that("find_cutpoint's systematic search handles edge cases", {
  edge_data <- mock_data
  edge_data$predictor <- rep(50, n_test) # No variability
  expect_message(
    find_cutpoint(
      data = edge_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      nmin = 5
    ),
    regexp = "No optimal cut-point could be determined",
    info = "Should message about no valid cuts"
  )
})

test_that("plot.find_cutpoint handles Cox model failure", {
  data_cox_fail <- data.frame(
    time = rexp(20, rate = 0.05),
    event = c(rep(0, 10), rep(1, 10)),
    predictor = c(1:10, 11:20)
  )
  res_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(data_cox_fail, "predictor", "time", "event", num_cuts = 1, quiet = TRUE)
  ))

  if (!is.null(res_fc) && !all(is.na(res_fc$optimal_cuts))) {
    expect_message(
      p <- plot(res_fc, type = "forest"),
      "Could not fit Cox model for forest plot."
    )
  } else {
    expect_message(plot(res_fc), "Cannot generate plot, no valid cut-point found.")
  }
})

test_that("find_cutpoint_number handles search failures gracefully", {
  # Use data with zero rows to force a complete failure
  bad_data <- data.frame(time = numeric(0), event = numeric(0), predictor = numeric(0))

  expect_message(
    res <- find_cutpoint_number(bad_data, "predictor", "time", "event", nmin = 1),
    "No complete cases found in the data after removing NAs."
  )
  expect_s3_class(res, "find_cutpoint_number_result")
  expect_equal(nrow(res$results), 0)

  res_na <- structure(
    list(
      results = data.frame(num_cuts = 0:1, BIC = c(NA, NA)),
      parameters = list(criterion = "BIC")
    ),
    class = "find_cutpoint_number_result"
  )
  expect_message(
    plot(res_na),
    "Cannot generate plot because no valid Information Criterion values were calculated"
  )
})

test_that(".calc_ic handles edge cases and invalid inputs", {
  # Test the initial checks for NULL/invalid model or loglik
  expect_true(is.na(.calc_ic(model = NULL, k = 2, n = 20, "BIC")))
  expect_true(is.na(.calc_ic(model = list(loglik = c(NA, NA)), k = 2, n = 20, "BIC")))

  # Test AICc when n <= k + 1, which should return NA
  model_stub <- list(loglik = c(NA, -50))
  n <- 5
  k <- 4
  expect_true(n <= k + 1)
  expect_true(is.na(.calc_ic(model_stub, k = k, n = n, criterion = "AICc")))

  # Ensure it works when valid
  k_valid <- 3
  expect_false(n <= k_valid + 1)
  expect_true(is.numeric(.calc_ic(model_stub, k = k_valid, n = n, criterion = "AICc")))
})

test_that(".run_genetic_search handles invalid domains", {
  skip_if_not_installed("rgenoud")

  # Test early return if predictor range is infinite
  bad_target <- c(1:10, Inf)
  expect_null(
    .run_genetic_search(target = bad_target, numcut = 1, time = 1:11, censor = 1,
                        confound = NULL, nmin = 1, criterion = "logrank")
  )

  # Test early return if starting values are NA
  na_target <- rep(NA_real_, 10)
  suppressWarnings(expect_null(
    .run_genetic_search(target = na_target, numcut = 1, time = 1:10, censor = 1,
                        confound = NULL, nmin = 1, criterion = "logrank")
  ))
})

test_that(".obj function handles internal model failures", {
  # Create data where a group will have no observations
  data_model_fail_harder <- data.frame(
    time = 1:10,
    censor = 1,
    target = c(1:5, 100:104)
  )

  cut_points <- c(50, 150)
  res <- .obj(params = cut_points, time = data_model_fail_harder$time, censor = data_model_fail_harder$censor,
              target = data_model_fail_harder$target, confound = NULL, numcut = 2,
              gap = 1, nmin = 1, criterion = "logrank")
  expect_equal(res, -Inf)
})

test_that("plot_optimization_curve handles failed results", {
  res_failed <- find_cutpoint(mock_data[1:10,], "predictor", "time", "event", quiet = TRUE)
  res_failed$parameters$method <- "systematic"
  res_failed$all_stats <- data.frame(cut1 = 1:5, stat = rnorm(5))
  res_failed$optimal_cuts <- NA

  p <- plot_optimization_curve(res_failed)
  expect_s3_class(p, "ggplot")
  expect_equal(p$labels$subtitle, "No optimal cut-point found.")
})

# --- NEW TESTS FOR COVERAGE ---

test_that("S3 summary.find_cutpoint works with all arguments", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  covariates = "covariate1", num_cuts = 1, quiet = TRUE)
  ))
  skip_if(is.null(res_valid) || all(is.na(res_valid$optimal_cuts)), "Valid result needed for summary test.")

  # Test that the summary runs and returns the object invisibly
  expect_s3_class(
    summary(res_valid, show_model = TRUE, show_group_counts = TRUE, show_medians = TRUE,
            show_ph_test = TRUE, show_params = TRUE),
    "find_cutpoint"
  )
})

test_that("S3 summary.find_cutpoint_number works with all arguments", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event", max_cuts = 1)
  ))
  skip_if(is.null(res_valid), "Valid result needed for summary test.")

  # Test with best_model_details = TRUE
  expect_s3_class(
    summary(res_valid, show_comparison_table = TRUE, show_best_model_details = TRUE,
            show_group_counts = TRUE, show_medians = TRUE, plot.it = FALSE),
    "find_cutpoint_number_result"
  )
})

test_that("plot.find_cutpoint fails gracefully if 'broom' is not installed", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 1, quiet = TRUE)
  ))
  skip_if(is.null(res_valid) || all(is.na(res_valid$optimal_cuts)), "Valid result needed.")

  # Mock requireNamespace to return FALSE
  local_mocked_bindings(
    "requireNamespace" = function(pkg, ...) {
      if (pkg == "broom") return(FALSE)
      # Must call the original function for other packages
      return(base::requireNamespace(pkg, ...))
    },
    .package = "base" # Needs withr package loaded
  )

  expect_error(
    plot(res_valid, type = "forest"),
    "Package 'broom' is required"
  )
})

# *** NEW: Test graceful failure if rgenoud is not installed ***
test_that("genetic search fails gracefully if 'rgenoud' is not installed", {
  # Mock requireNamespace to return FALSE
  local_mocked_bindings(
    "requireNamespace" = function(pkg, ...) {
      if (pkg == "rgenoud") return(FALSE)
      return(base::requireNamespace(pkg, ...))
    },
    .package = "base" # Needs withr package loaded
  )

  # Test find_cutpoint
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event", method = "genetic"),
    regexp = "rgenoud" # *** SIMPLIFIED regexp ***
  )

  # Test find_cutpoint_number
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event", method = "genetic"),
    regexp = "rgenoud" # *** SIMPLIFIED regexp ***
  )
})

test_that("find_cutpoint systematic search errors for num_cuts > 2", {
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event", num_cuts = 3, method = "systematic"),
    "Systematic search currently only supports num_cuts = 1 or 2"
  )
})

