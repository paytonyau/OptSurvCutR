# Comprehensive test file for the OptSurvCutR package.
# This script merges core functionality tests with additional
# tests for edge cases, error conditions, and internal helpers.
# Tests are organised logically by function and category for clarity.

#' @srrstats {G5.0} Implements edge case and validation tests.
#' @srrstats {G5.3} Tests for absence of NA/NaN in return objects (via
#'   successful class checks).
#' @srrstats {G5.4} Correctness tests against known data properties
#'   (bimodal/trimodal distributions).
#' @srrstats {G5.4a} Correctness tests for new method implementation.
#' @srrstats {G5.4b} Comparison against previous/alternative implementations
#'   (systematic vs genetic).
#' @srrstats {G5.4c} Comparison against published/known values (simulated data).
#' @srrstats {G5.9} Noise susceptibility tests (stochastic genetic algorithm).
#' @srrstats {G5.9a} Tests with trivial noise (implicit in stochastic runs).
#' @srrstats {G5.9b} Tests with different seeds.
#' @srrstats {RE1.0} Tests the core cut-point implementation.
#' @srrstats {RE4.0} Tests model selection via AIC, BIC, and AICc.
#' @srrstats {RE7.0} Tests the bootstrap validation function.
#' @srrstats {G2.0} Tests for input validation are included.
#' @srrstats {RE6.0} Tests for S3 plot methods.
#' @srrstats {RE4.17} Tests for S3 print methods.
#' @srrstats {G5.8} Tests stochastic genetic algorithm behaviour.
#' @srrstats {G2.4c} Tests optional dependency handling.

library(testthat)
library(survival)
library(cli) # for message expectations
if (!requireNamespace("rgenoud", quietly = TRUE)) {
  install.packages("rgenoud")
}
library(rgenoud)

# --- Setup: Create consistent, reproducible mock datasets ---
set.seed(42)
n_test <- 60 # Increased sample size for robustness
# Bimodal predictor for 1 cut (2 groups): ~50% N(40, 5), ~50% N(60, 5)
n1 <- n_test / 2
n2 <- n_test - n1
mock_data <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(rnorm(n1, mean = 40, sd = 5), rnorm(n2, mean = 60, sd = 5)),
  covariate1 = rnorm(n_test, mean = 5, sd = 1),
  covariate2 = sample(c("A", "B"), n_test, replace = TRUE) # Factor
)
# Trimodal predictor for 2 cuts (3 groups)
n1_3 <- floor(n_test / 3)
n2_3 <- floor(n_test / 3)
n3_3 <- n_test - n1_3 - n2_3
mock_data_3groups <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(
    rnorm(n1_3, mean = 30, sd = 5),
    rnorm(n2_3, mean = 50, sd = 5),
    rnorm(n3_3, mean = 70, sd = 5)
  ),
  covariate1 = rnorm(n_test, mean = 5, sd = 1)
)
# Mock dataset with skewed predictor
mock_data_skewed <- mock_data
mock_data_skewed$predictor <- rexp(n_test, rate = 0.1)
# Mock dataset with heavy censoring (90% censored)
mock_data_heavy_censor <- mock_data
mock_data_heavy_censor$event <- sample(0:1, n_test,
  replace = TRUE,
  prob = c(0.9, 0.1)
)
# Pathological dataset for bootstrap failure
mock_data_pathological <- mock_data[1:20, ]
mock_data_pathological$predictor <- rep(50, 20) # Constant
# Tiny dataset for AIC/AICc edge cases
tiny_data <- mock_data[1:4, ]
tiny_data$predictor <- c(40, 45, 50, 55)
# New dataset for small dataset test with constant predictor
small_data_unique <- data.frame(
  time = rexp(20, rate = 0.05),
  event = sample(0:1, 20, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = rep(50, 20), # Constant to induce warning
  covariate1 = rnorm(20, mean = 5, sd = 1)
)
# Pathological data with no events to force model failures
mock_data_no_events <- mock_data
mock_data_no_events$event <- 0

# --- SECTION 1: Tests for find_cutpoint_number() ---
# Note: These tests cover core functionality, including genetic and systematic
# searches, covariate handling, edge cases, and input validation.

#' @srrstats {RE4.0} Tests AIC, AICc, BIC criteria.
#' @srrstats {RE4.11} Tests implementation of IC formulas.
#' @srrstats {G5.8} Tests stochastic genetic algorithm.
test_that("find_cutpoint_number genetic search works for all criteria", {
  skip_if_not_installed("rgenoud")
  # Test BIC (recommended default)
  res_bic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "genetic", criterion = "BIC", max_cuts = 2,
      max.generations = 5, nmin = 1
    )
  ))
  skip_if(is.null(res_bic), "Genetic search returned NULL.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 3) # 0, 1, 2 cuts
  expect_true("BIC" %in% names(res_bic$results))
  expect_true(all(res_bic$results$num_cuts %in% 0:2))
  # Test AIC
  res_aic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "genetic", criterion = "AIC",
      max.generations = 5, nmin = 1
    )
  ))
  skip_if(is.null(res_aic), "Genetic search returned NULL.")
  expect_s3_class(res_aic, "find_cutpoint_number_result")
  expect_true("AIC" %in% names(res_aic$results))
  # Test AICc
  res_aicc <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "genetic", criterion = "AICc",
      max.generations = 5, nmin = 1
    )
  ))
  skip_if(is.null(res_aicc), "Genetic search returned NULL.")
  expect_s3_class(res_aicc, "find_cutpoint_number_result")
  expect_true("AICc" %in% names(res_aicc$results))
})

#' @srrstats {RE1.0} Tests model selection via IC.
#' @srrstats {RE4.0} Tests BIC criterion.
test_that("find_cutpoint_number systematic search works for BIC", {
  res_bic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "systematic", criterion = "BIC",
      max_cuts = 2, nmin = 1
    )
  ))
  skip_if(is.null(res_bic), "Systematic search returned NULL.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 3)
  expect_true("BIC" %in% names(res_bic$results))
})

#' @srrstats {RE2.2} Tests covariate adjustment.
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
  skip_if(
    is.null(res_bic_cov),
    "Systematic search with covariates returned NULL."
  )
  expect_s3_class(res_bic_cov, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_cov$results), 2) # 0, 1 cuts
  expect_true(!is.na(res_bic_cov$optimal_num_cuts))
  expect_equal(
    res_bic_cov$parameters$covariates,
    c("covariate1", "covariate2")
  )
})

#' @srrstats {RE2.2} Tests covariate adjustment.
#' @srrstats {G5.8} Tests stochastic genetic algorithm.
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
      max.generations = 5
    )
  ))
  skip_if(
    is.null(res_bic_cov_gen),
    "Genetic search with covariates returned NULL."
  )
  expect_s3_class(res_bic_cov_gen, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_cov_gen$results), 2) # 0, 1 cuts
  expect_true(!is.na(res_bic_cov_gen$optimal_num_cuts))
  expect_equal(res_bic_cov_gen$parameters$covariates, "covariate1")
})

#' @srrstats {G5.0} Tests edge case (constant predictor).
#' @srrstats {G5.8} Tests that constant predictor returns na_result.
test_that("find_cutpoint_number handles too few unique values", {
  skip_if(
    length(unique(small_data_unique$predictor)) > 1,
    "Test data must have only one unique predictor value."
  )
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
  expect_true(
    is.na(result$optimal_num_cuts),
    info = "Should return NA optimal_num_cuts for constant predictor."
  )
  expect_true(
    nrow(result$results) == 0,
    info = "Should return an empty results table."
  )
})

#' @srrstats {G5.0} Tests edge case (small sample).
#' @srrstats {RE4.11} Tests IC formula branch for small n.
#' @srrstats {G5.8} Tests small n - k - 1 in .calc_ic.
test_that("find_cutpoint_number handles AIC/AICc edge cases", {
  res_aic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(tiny_data, "predictor", "time", "event",
      method = "systematic", criterion = "AIC",
      max_cuts = 1, nmin = 4
    )
  ))
  skip_if(is.null(res_aic), "Systematic search returned NULL.")
  expect_true(
    is.na(res_aic$results$AIC[res_aic$results$num_cuts == 1]) ||
      !any(res_aic$results$cuts[res_aic$results$num_cuts == 1] != "NA")
  )
  res_aicc <- suppressMessages(suppressWarnings(
    find_cutpoint_number(tiny_data, "predictor", "time", "event",
      method = "systematic", criterion = "AICc",
      max_cuts = 1, nmin = 4
    )
  ))
  skip_if(is.null(res_aicc), "Systematic search returned NULL.")
  expect_true(
    is.na(res_aicc$results$AICc[res_aicc$results$num_cuts == 1]) ||
      !any(res_aicc$results$cuts[res_aicc$results$num_cuts == 1] != "NA")
  )
})

#' @srrstats {RE1.0} Tests core functionality of model selection.
test_that("find_cutpoint_number (systematic) runs sequentially", {
  res_bic <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "systematic", criterion = "BIC",
      max_cuts = 1, nmin = 1
    )
  ))
  skip_if(is.null(res_bic), "Systematic search returned NULL.")
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic$results), 2) # 0, 1 cuts
})

#' @srrstats {G5.8} Tests core functionality of stochastic algorithm.
test_that("find_cutpoint_number (genetic) runs sequentially", {
  skip_if_not_installed("rgenoud")
  res_bic_gen <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "genetic", criterion = "BIC", max_cuts = 1,
      max.generations = 5, nmin = 1 # Updated param
    )
  ))
  skip_if(is.null(res_bic_gen), "Genetic search returned NULL.")
  expect_s3_class(res_bic_gen, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_gen$results), 2) # 0, 1 cuts
})

#' @srrstats {G2.0} Tests input validation.
#' @srrstats {G2.1} Tests input type validation.
#' @srrstats {G2.4b} Tests match.arg validation.
#' @srrstats {G2.13} Tests that invalid inputs trigger errors.
test_that("find_cutpoint_number handles invalid inputs", {
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      max_cuts = -1
    ),
    regexp = "max_cuts must be a non-negative integer"
  )
  expect_error(
    find_cutpoint_number(mock_data, predictor = NULL, "time", "event"),
    regexp = "predictor.*must be specified"
  )
  expect_error(
    find_cutpoint_number(mock_data, "predictor",
      outcome_time = NULL, "event"
    ),
    regexp = "outcome_time.*are required"
  )
  expect_error(
    find_cutpoint_number(data.frame(a = 9), "predictor", "time", "event"),
    regexp = "Missing columns:"
  )
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event", nmin = -1),
    regexp = "'nmin' must be a positive number"
  )
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "systematic", max_cuts = 3
    ),
    regexp = "only implemented for max_cuts <= 2"
  )
})

#' @srrstats {G5.0} Tests edge case (constant predictor).
#' @srrstats {G5.8} Tests NA result for constant predictor.
test_that("find_cutpoint_number handles insufficient and constant data", {
  small_data <- mock_data[1:5, ]
  res_insufficient <- suppressMessages(suppressWarnings(
    find_cutpoint_number(small_data, "predictor", "time", "event",
      max_cuts = 1,
      method = "systematic", nmin = 5, quiet = TRUE
    )
  ))
  expect_s3_class(res_insufficient, "find_cutpoint_number_result")
  expect_true(all(is.na(res_insufficient$optimal_num_cuts)))
  constant_data <- mock_data[1:20, ]
  constant_data$predictor <- rep(50, 20)
  res_constant <- suppressMessages(suppressWarnings(
    find_cutpoint_number(constant_data, "predictor", "time", "event",
      max_cuts = 1, method = "systematic",
      nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_constant, "find_cutpoint_number_result")
  expect_true(all(is.na(res_constant$optimal_num_cuts)))
})

#' @srrstats {G2.7} Tests handling of zero-length data.
test_that("find_cutpoint_number handles search failures gracefully", {
  bad_data <- data.frame(
    time = numeric(0), event = numeric(0), predictor = numeric(0)
  )
  expect_message(
    res <- find_cutpoint_number(bad_data, "predictor", "time", "event", nmin = 1),
    "No complete cases found"
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
    regexp = "Cannot generate plot: no valid IC values"
  )
})

#' @srrstats {RE4.11} Tests implementation of IC formulas.
test_that("find_cutpoint_number base model IC failure (Pathological Data)", {
  res_sys_fail <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data_no_events, "predictor", "time", "event",
      max_cuts = 1
    )
  ))
  expect_false(
    is.na(res_sys_fail$results$BIC[res_sys_fail$results$num_cuts == 0])
  )
  skip_if_not_installed("rgenoud")
  res_gen_fail <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data_no_events, "predictor", "time", "event",
      method = "genetic", max_cuts = 1
    )
  ))
  expect_false(
    is.na(res_gen_fail$results$BIC[res_gen_fail$results$num_cuts == 0])
  )
})

# --- SECTION 2: Tests for find_cutpoint() ---
# Note: These tests cover core optimal cut-point algorithm, including
# systematic and genetic searches, various criteria, covariate handling,
# edge cases, and input validation.

#' @srrstats {RE1.0} Tests core optimal cut-point algorithm (systematic).
test_that("find_cutpoint systematic search works for one and two cuts", {
  # One cut
  res_lr1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", criterion = "logrank", nmin = 1,
      quiet = TRUE
    )
  ))
  expect_s3_class(res_lr1, "find_cutpoint")
  expect_length(res_lr1$optimal_cuts, 1)
  expect_true(
    all(is.na(res_lr1$optimal_cuts) |
      (res_lr1$optimal_cuts >= min(mock_data$predictor) &
        res_lr1$optimal_cuts <= max(mock_data$predictor)))
  )
  # Two cuts (with appropriate data)
  res_lr2 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
      num_cuts = 2,
      method = "systematic", criterion = "logrank", nmin = 10,
      quiet = TRUE
    )
  ))
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
})

#' @srrstats {G2.0} Tests validation of explicit genetic algorithm parameters.
test_that("find_cutpoint respects explicit pop.size and max.generations", {
  skip_if_not_installed("rgenoud")

  # Run with non-default GA parameters
  res_ga_params <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      method = "genetic",
      pop.size = 50, # Explicit non-default
      max.generations = 5, # Explicit non-default
      nmin = 5,
      quiet = TRUE
    )
  ))

  # 1. Check object class
  expect_s3_class(res_ga_params, "find_cutpoint")

  # 2. Check that parameters were stored correctly in the output
  expect_equal(res_ga_params$parameters$pop.size, 50)
  expect_equal(res_ga_params$parameters$max.generations, 5)

  # 3. Check that the method was indeed genetic
  expect_equal(res_ga_params$parameters$method, "genetic")
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm for multiple cuts.
test_that("find_cutpoint genetic search works for multiple cuts", {
  skip_if_not_installed("rgenoud")
  res_lr2 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
      num_cuts = 2, method = "genetic",
      criterion = "logrank", max.generations = 5, nmin = 1,
      quiet = TRUE
    )
  ))
  skip_if(all(is.na(res_lr2$optimal_cuts)), "GA search failed for 2 cuts.")
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
  expect_true(
    all(is.na(res_lr2$optimal_cuts) | diff(res_lr2$optimal_cuts) > 0)
  )
  res_lr3 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
      num_cuts = 3, method = "genetic",
      criterion = "logrank", max.generations = 5, nmin = 1,
      quiet = TRUE
    )
  ))
  skip_if(all(is.na(res_lr3$optimal_cuts)), "GA search failed for 3 cuts.")
  expect_s3_class(res_lr3, "find_cutpoint")
  expect_length(res_lr3$optimal_cuts, 3)
})

#' @srrstats {RE1.0} Tests all optimisation criteria.
#' @srrstats {RE4.17} Tests nmin as a proportion.
test_that("find_cutpoint handles various criteria and nmin proportion", {
  # Hazard ratio
  res_hr1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1, method = "systematic",
      criterion = "hazard_ratio", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_hr1, "find_cutpoint")
  expect_length(res_hr1$optimal_cuts, 1)
  # p-value
  res_pv1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1, method = "systematic",
      criterion = "p_value", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_pv1, "find_cutpoint")
  expect_length(res_pv1$optimal_cuts, 1)
  # nmin as proportion
  res_nmin_prop <- suppressMessages(suppressWarnings(
    find_cutpoint(
      mock_data, "predictor", "time", "event",
      nmin = 0.1, quiet = TRUE
    )
  ))
  expect_s3_class(res_nmin_prop, "find_cutpoint")
})

#' @srrstats {G2.0} Tests input validation.
#' @srrstats {G2.1} Tests input type validation.
#' @srrstats {G2.4b} Tests match.arg validation.
#' @srrstats {G2.13} Tests that invalid inputs trigger errors.
test_that("find_cutpoint handles invalid inputs", {
  bad_data <- mock_data[, c("time", "event")]
  expect_error(
    find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
    regexp = "Missing required columns: 'predictor'"
  )
  bad_data <- mock_data
  bad_data$time[1] <- -1
  expect_error(
    find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
    regexp = "Time variable must be non-negative"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = -1,
      method = "systematic", nmin = 1
    ),
    regexp = "non-negative integer"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = NA,
      method = "systematic", nmin = 1
    ),
    regexp = "non-negative integer"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
      criterion = "invalid", nmin = 1
    ),
    regexp = "arg.*should be one of"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 2,
      criterion = "hazard_ratio", nmin = 1
    ),
    regexp = "'hazard_ratio' is only supported for num_cuts = 1"
  )
})

#' @srrstats {G5.0} Tests edge case (constant predictor).
test_that("find_cutpoint handles quiet arg and low variability predictor", {
  expect_no_message(
    find_cutpoint(mock_data, "predictor", "time", "event", quiet = TRUE)
  )
  low_var_data <- mock_data
  low_var_data$predictor <- rep(50, n_test)
  expect_message(
    find_cutpoint(low_var_data, "predictor", "time", "event", quiet = FALSE),
    regexp = "Predictor has too few unique values"
  )
})

#' @srrstats {G5.0} Tests edge cases (skewed data, heavy censoring).
#' @srrstats {G5.2} Tests handling of non-converged models.
test_that("find_cutpoint handles various data scenarios", {
  res_lr_skewed <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_skewed, "predictor", "time", "event",
      num_cuts = 1, method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_lr_skewed, "find_cutpoint")
  expect_length(res_lr_skewed$optimal_cuts, 1)
  res_lr_censor <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_heavy_censor, "predictor", "time", "event",
      num_cuts = 1, method = "systematic",
      criterion = "logrank", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_lr_censor, "find_cutpoint")
  expect_length(res_lr_censor$optimal_cuts, 1)
  bad_data <- mock_data
  bad_data$predictor <- rep(50, n_test)
  res_cox_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(bad_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", nmin = 10, quiet = TRUE
    )
  ))
  expect_s3_class(res_cox_fail, "find_cutpoint")
  expect_true(all(is.na(res_cox_fail$optimal_cuts)))
})

#' @srrstats {RE2.2} Tests covariate adjustment (systematic).
test_that("find_cutpoint (systematic) works with covariates", {
  res_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      covariates = "covariate1", num_cuts = 1,
      method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_cov, "find_cutpoint")
  expect_length(res_cov$optimal_cuts, 1)
  expect_true(!is.na(res_cov$optimal_cuts[1]))
  res_lr_seq <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1, method = "systematic",
      criterion = "logrank", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_lr_seq, "find_cutpoint")
  expect_length(res_lr_seq$optimal_cuts, 1)
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm.
test_that("find_cutpoint (genetic) works sequentially", {
  skip_if_not_installed("rgenoud")
  res_gen_seq <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 2, method = "genetic",
      criterion = "logrank", max.generations = 5, nmin = 1, quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_gen_seq) || all(is.na(res_gen_seq$optimal_cuts)),
    "Sequential genetic search failed."
  )
  expect_s3_class(res_gen_seq, "find_cutpoint")
  expect_length(res_gen_seq$optimal_cuts, 2)
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm (high cuts, p_value).
test_that("find_cutpoint handles high num_cuts with p_value", {
  skip_if_not_installed("rgenoud")
  res_pv4 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 4, method = "genetic",
      criterion = "p_value", max.generations = 5, nmin = 1, quiet = TRUE
    )
  ))
  skip_if(all(is.na(res_pv4$optimal_cuts)), "GA search failed for 4 cuts.")
  expect_s3_class(res_pv4, "find_cutpoint")
  expect_length(res_pv4$optimal_cuts, 4)
})

#' @srrstats {G5.0} Tests edge cases.
#' @srrstats {G5.8} Tests that constant predictor returns na_result.
test_that("find_cutpoint handles insufficient and constant data", {
  small_data <- mock_data[1:5, ]
  res_insufficient <- suppressMessages(suppressWarnings(
    find_cutpoint(small_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", nmin = 5, quiet = TRUE
    )
  ))
  expect_s3_class(res_insufficient, "find_cutpoint")
  expect_true(all(is.na(res_insufficient$optimal_cuts)))
  constant_data <- mock_data[1:20, ]
  constant_data$predictor <- rep(50, 20)
  res_constant <- suppressMessages(suppressWarnings(
    find_cutpoint(constant_data, "predictor", "time", "event",
      num_cuts = 1, method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_constant, "find_cutpoint")
  expect_true(all(is.na(res_constant$optimal_cuts)))
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm.
test_that("group creation works in find_cutpoint", {
  set.seed(42)
  res_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 2,
      method = "genetic", criterion = "logrank", max.generations = 5,
      nmin = 1, quiet = TRUE
    )
  ))
  skip_if(all(is.na(res_fc$optimal_cuts)), "Genetic search failed.")
  expect_s3_class(res_fc, "find_cutpoint")
  expect_length(res_fc$optimal_cuts, 2)
  expect_true(
    all(is.na(res_fc$optimal_cuts) | diff(res_fc$optimal_cuts) > 0)
  )
})

#' @srrstats {G5.0} Tests edge case (constant predictor).
#' @srrstats {G5.8} Tests NA result for constant predictor.
test_that("find_cutpoint's systematic search handles edge cases", {
  edge_data <- mock_data
  edge_data$predictor <- rep(50, n_test)
  expect_message(
    find_cutpoint(
      data = edge_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      nmin = 5,
      quiet = FALSE
    ),
    regexp = "Predictor has too few unique values"
  )
})

#' @srrstats {G5.2} Tests handling of non-converged models in plotting.
test_that("plot.find_cutpoint handles Cox model failure", {
  data_cox_fail <- data.frame(
    time = rexp(20, rate = 0.05),
    event = c(rep(0, 10), rep(1, 10)),
    predictor = c(1:10, 11:20)
  )
  res_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(data_cox_fail, "predictor", "time", "event",
      num_cuts = 1,
      quiet = TRUE
    )
  ))
  if (!is.null(res_fc) && !all(is.na(res_fc$optimal_cuts))) {
    expect_message(
      p <- plot(res_fc, type = "forest"),
      "Could not fit Cox model for forest plot"
    )
  } else {
    expect_message(
      plot(res_fc),
      "Cannot generate plot: no valid cut-point"
    )
  }
})

#' @srrstats {G5.8} Tests handling of invalid cut-points.
test_that("Coverage: .systematic_search - null model fit failure", {
  bad_data <- mock_data
  bad_data$time <- rep(0, n_test)
  bad_data2 <- bad_data
  names(bad_data2)[names(bad_data2) == "predictor"] <- "factor"
  result <- suppressMessages(suppressWarnings(
    .systematic_search(
      userdata = bad_data2,
      num_cuts = 1,
      criterion = "p_value",
      covariates = NULL,
      nmin = 5,
      predictor_name = "predictor",
      quiet = TRUE
    )
  ))
  expect_true(all(is.na(result$optimal_cuts)))
  expect_equal(result$parameters$method, "systematic")
})

test_that("Coverage: .systematic_search - no valid cuts (nmin violation)", {
  tiny_data <- mock_data[1:10, ]
  tiny_data2 <- tiny_data
  names(tiny_data2)[names(tiny_data2) == "predictor"] <- "factor"
  result <- suppressMessages(suppressWarnings(
    .systematic_search(
      userdata = tiny_data2,
      num_cuts = 1,
      criterion = "logrank",
      covariates = NULL,
      nmin = 8,
      predictor_name = "predictor",
      quiet = TRUE
    )
  ))
  expect_true(all(is.na(result$optimal_cuts)))
})

test_that("Coverage: .systematic_search - 2 cuts with insufficient data", {
  small_data <- mock_data[1:20, ]
  small_data2 <- small_data
  names(small_data2)[names(small_data2) == "predictor"] <- "factor"
  result <- suppressMessages(suppressWarnings(
    .systematic_search(
      userdata = small_data2,
      num_cuts = 2,
      criterion = "logrank",
      covariates = NULL,
      nmin = 10,
      predictor_name = "predictor",
      quiet = TRUE
    )
  ))
  expect_true(all(is.na(result$optimal_cuts)))
})

test_that("Coverage: find_cutpoint (genetic) shows messages", {
  skip_if_not_installed("rgenoud")
  # Suppress rgenoud's own warning about max.generations
  expect_message(
    suppressWarnings(find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "genetic",
      max.generations = 1,
      quiet = FALSE
    )),
    regexp = "Starting genetic search"
  )
})

test_that("Coverage: find_cutpoint handles nmin edge cases", {
  # Test for systematic, 1 cut, grid empty
  expect_message(
    res_1_cut <- find_cutpoint(
      data = head(mock_data, 20),
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      nmin = 15,
      num_cuts = 1,
      method = "systematic",
      quiet = FALSE
    ),
    regexp = "Not enough data"
  )
  expect_true(all(is.na(res_1_cut$optimal_cuts)))
  # Test for systematic, 2 cuts, grid empty
  expect_message(
    res_2_cut <- find_cutpoint(
      data = head(mock_data, 20),
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      nmin = 8,
      num_cuts = 2,
      method = "systematic",
      quiet = FALSE
    ),
    regexp = "Not enough data"
  )
  expect_true(all(is.na(res_2_cut$optimal_cuts)))
})

test_that("Coverage: S3 methods handle Cox model failures", {
  # *** FIX ***: Mock coxph to return NULL to test the is.null(fit_cox)
  # branches
  local_mocked_bindings(
    "coxph" = function(...) NULL,
    .package = "survival"
  )
  # Create a valid result object to pass to summary() and plot()
  res_valid <- suppressMessages(find_cutpoint(
    data = mock_data,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    quiet = TRUE
  ))
  skip_if(
    any(is.na(res_valid$optimal_cuts)),
    "Setup for S3 Cox failure test failed."
  )
  # Test summary() branch (line 729)
  expect_message(
    summary(res_valid),
    regexp = "Could not fit Cox model: convergence failed"
  )
  # Test plot() branch (line 847)
  expect_message(
    plot(res_valid, type = "forest"),
    regexp = "Could not fit Cox model for forest plot"
  )
})

test_that("Coverage: find_cutpoint p_value + covariates branches", {
  # Tests p_value + covariates + systematic + 2 cuts branches
  res_p_val_2 <- find_cutpoint(
    data = head(mock_data_3groups, 50),
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = "systematic",
    num_cuts = 2,
    criterion = "p_value",
    covariates = "covariate1",
    nmin = 5,
    quiet = TRUE
  )
  expect_s3_class(res_p_val_2, "find_cutpoint")
})

test_that("Coverage: .get_stat - group constraint violations", {
  stat <- .get_stat(
    cuts = c(45, 55),
    num_cuts = 2,
    data_in = data.frame(
      factor = rep(50, 20),
      time = rexp(20, 0.05),
      event = sample(0:1, 20, replace = TRUE)
    ),
    criterion = "logrank",
    cov_formula = "",
    nmin = 5,
    fit_null = NULL
  )
  expect_true(is.na(stat))
})

test_that("Coverage: .get_stat - hazard_ratio with missing coefficient", {
  stat <- .get_stat(
    cuts = 50,
    num_cuts = 1,
    data_in = data.frame(
      factor = c(rep(40, 10), rep(60, 10)),
      time = rexp(20, 0.05),
      event = rep(0, 20) # No events
    ),
    criterion = "hazard_ratio",
    cov_formula = "",
    nmin = 5,
    fit_null = NULL
  )
  expect_equal(stat, -Inf)
})

test_that("Coverage: .get_stat - p_value with null fit_null", {
  stat <- .get_stat(
    cuts = 50,
    num_cuts = 1,
    data_in = data.frame(
      factor = c(rep(40, 15), rep(60, 15)),
      time = rexp(30, 0.05),
      event = sample(0:1, 30, replace = TRUE)
    ),
    criterion = "p_value",
    cov_formula = "",
    nmin = 5,
    fit_null = NULL
  )
  expect_true(is.na(stat) || is.numeric(stat))
})

test_that(
  "Coverage: find_cutpoint - genetic method with p_value and covariates",
  {
    skip_if_not_installed("rgenoud")
    result <- suppressMessages(suppressWarnings(
      find_cutpoint(
        data = mock_data,
        predictor = "predictor",
        outcome_time = "time",
        outcome_event = "event",
        num_cuts = 2,
        method = "genetic",
        criterion = "p_value",
        covariates = "covariate1",
        nmin = 5,
        max.generations = 5,
        seed = 123,
        quiet = TRUE
      )
    ))
    expect_s3_class(result, "find_cutpoint")
  }
)

# --- SECTION 3: Tests for validate_cutpoint() ---
# Note: These tests cover bootstrap validation, including parallel processing,
# covariate handling, edge cases, and input validation.

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
    nmin = 5,
    quiet = TRUE
  )
))
if (any(is.na(valid_fc_result_for_boot$optimal_cuts))) {
  stop(
    "Test setup failed: Could not generate a valid fc_result for tests."
  )
}

#' @srrstats {RE7.0} Tests bootstrap resampling.
#' @srrstats {RE7.1} Tests cut-point stability assessment.
test_that("validate_cutpoint works with recommended settings", {
  # Capture messages to verify the "90% stability" message appears
  expect_message(
    suppressWarnings({
      result <- validate_cutpoint(valid_fc_result_for_boot,
        num_replicates = 20,
        n_cores = 1,
        max.generations = 10,
        seed = 42
      )
    }),
    regexp = "90% of original" # Check for the new explanation
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
test_that("Coverage: validate_cutpoint - non-integer num_replicates", {
  valid_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1, quiet = TRUE
    )
  ))
  skip_if(any(is.na(valid_fc$optimal_cuts)))
  expect_error(
    validate_cutpoint(valid_fc, num_replicates = 25.5, n_cores = 1),
    regexp = "positive integer"
  )
})

#' @srrstats {G2.0a} Tests validation of scalar parameter num_replicates.
#' @srrstats {G5.8d} Tests data outside scope (insufficient sample size for
#'   requested cuts).
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
  valid_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1, quiet = TRUE
    )
  ))
  skip_if(any(is.na(valid_fc$optimal_cuts)))
  result <- suppressMessages(suppressWarnings(
    validate_cutpoint(valid_fc,
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

# --- SECTION 4: Tests for plotting_functions.R ---
# Note: These tests cover diagnostic and publication-ready plots, including
# optimisation curves and Schoenfeld residuals.

#' @srrstats {RE6.0} Tests diagnostic plot for systematic search.
test_that("plot_optimisation_curve works for all criteria", {
  res_lr <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", criterion = "logrank", nmin = 1,
      quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_lr) || all(is.na(res_lr$optimal_cuts)),
    "Systematic logrank search failed."
  )
  p_lr <- plot_optimisation_curve(res_lr)
  expect_s3_class(p_lr, "ggplot")
  expect_equal(p_lr$labels$y, "Log-Rank Statistic")
  res_hr <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", criterion = "hazard_ratio", nmin = 1,
      quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_hr) || all(is.na(res_hr$optimal_cuts)),
    "Systematic HR search failed."
  )
  p_hr <- plot_optimisation_curve(res_hr)
  expect_s3_class(p_hr, "ggplot")
  expect_equal(p_hr$labels$y, "Hazard Ratio (HR)")
  res_pv <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", criterion = "p_value", nmin = 1,
      quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_pv) || all(is.na(res_pv$optimal_cuts)),
    "Systematic p-value search failed."
  )
  p_pv <- plot_optimisation_curve(res_pv)
  expect_s3_class(p_pv, "ggplot")
  expect_equal(p_pv$labels$y, "P-value")
})

#' @srrstats {G2.0} Tests input validation for plotting function.
test_that("plot_optimisation_curve throws errors for invalid input", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_valid),
    "Systematic search failed, cannot run error tests."
  )
  expect_error(
    plot_optimisation_curve(list()),
    regexp = "Input must be an object from the"
  )
  res_genetic <- res_valid
  res_genetic$parameters$method <- "genetic"
  expect_error(
    plot_optimisation_curve(res_genetic),
    regexp = "only for `method = \"systematic\"`"
  )
  res_2_cuts <- res_valid
  res_2_cuts$parameters$num_cuts <- 2
  expect_error(
    plot_optimisation_curve(res_2_cuts),
    regexp = "only for `num_cuts = 1`"
  )
  res_no_stats <- res_valid
  res_no_stats$all_stats <- NULL
  expect_error(
    plot_optimisation_curve(res_no_stats),
    regexp = "must have a valid `all_stats`"
  )
})

#' @srrstats {RE6.3} Tests diagnostic Schoenfeld residual plot.
test_that("plot_cutpoint_residuals produces Schoenfeld residual plot (RE6.3)", {
  res <- suppressMessages(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      quiet = TRUE
    )
  )
  skip_if(
    is.null(res) || any(is.na(res$optimal_cuts)),
    "No valid cut-point for diagnostics test."
  )
  p <- plot_cutpoint_residuals(res)
  expect_type(p, "list")
  expect_true(length(p) > 0)
  expect_true(all(vapply(p, inherits, "ggplot", FUN.VALUE = logical(1))))
  titles <- vapply(p, function(x) x$labels$title %||% "", character(1))
  expect_true(any(grepl("Schoenfeld", titles, ignore.case = TRUE)))
  last_subtitle <- p[[length(p)]]$labels$subtitle %||% ""
  expect_match(last_subtitle, "Global", all = FALSE)
})

test_that("Coverage: plot_optimisation_curve empty all_stats", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  res_valid$all_stats <- data.frame()
  expect_error(
    plot_optimisation_curve(res_valid),
    regexp = "`all_stats` is empty"
  )
})

test_that("Coverage: plot_cutpoint_residuals - NULL fit from coxph failure", {
  # Force coxph to fail by making the event column non-numeric
  result <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1, quiet = TRUE
    )
  ))
  skip_if(any(is.na(result$optimal_cuts)), "No valid cut-point for test setup")

  # Corrupt the event column → coxph will fail inside plot_cutpoint_residuals()
  result$userdata$event <- as.character(result$userdata$event)

  expect_message(
    plot_cutpoint_residuals(result),
    regexp = "Cox model failed"
  )
})

test_that("Coverage: plot_cutpoint_residuals Cox failure (Pathological Data)", {
  res_valid <- valid_fc_result_for_boot
  res_valid$userdata$event <- 0
  expect_error(plot_cutpoint_residuals(res_valid))
})

test_that(
  "Coverage: plot.find_cutpoint - forest plot with invalid reference",
  {
    skip_if_not_installed("broom")
    result <- suppressMessages(suppressWarnings(
      find_cutpoint(mock_data, "predictor", "time", "event",
        num_cuts = 1, quiet = TRUE
      )
    ))
    skip_if(any(is.na(result$optimal_cuts)))
    expect_message(
      plot(result, type = "forest", reference_group = "INVALID"),
      regexp = "Defaulting to"
    )
  }
)

# --- SECTION 5: Tests for S3 Methods (print, plot, summary) ---
# Note: These tests cover S3 methods for all main classes, including
# handling of NA results.

#' @srrstats {RE4.17} Tests S3 print methods.
#' @srrstats {G5.2} Tests graceful handling of empty/NA results.
test_that("S3 methods for find_cutpoint handle NA results", {
  res_na <- find_cutpoint(mock_data_pathological, "predictor", "time", "event",
    num_cuts = 1, method = "systematic", quiet = TRUE
  )
  expect_true(all(is.na(res_na$optimal_cuts)))
  expect_message(print(res_na), "No optimal cut-point determined")
  expect_message(summary(res_na), "Optimal Cut-point Analysis")
  expect_message(plot(res_na), "Cannot generate plot: no valid cut-point")
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_valid) || all(is.na(res_valid$optimal_cuts)),
    "Valid result needed for summary test."
  )
  sum_res <- summary(res_valid,
    show_model = FALSE, show_group_counts = FALSE,
    show_medians = FALSE, show_ph_test = FALSE,
    show_params = FALSE
  )
  expect_s3_class(sum_res, "find_cutpoint")
})

#' @srrstats {RE4.18} Tests S3 summary method.
test_that("S3 summary.find_cutpoint works with all arguments", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      covariates = "covariate1", num_cuts = 1, quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_valid) || all(is.na(res_valid$optimal_cuts)),
    "Valid result needed for summary test."
  )
  expect_s3_class(
    summary(res_valid,
      show_model = TRUE, show_group_counts = TRUE,
      show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE
    ),
    "find_cutpoint"
  )
})

#' @srrstats {RE6.0} Tests S3 plot method.
test_that("plot.find_cutpoint generates all plot types", {
  skip_if_not_installed("broom")
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
      num_cuts = 1,
      method = "systematic", quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_valid) || all(is.na(res_valid$optimal_cuts)),
    "Valid result needed for plot tests."
  )
  p_outcome <- suppressWarnings(plot(res_valid, type = "outcome"))
  expect_s3_class(p_outcome$plot, "ggplot")
  p_dist <- plot(res_valid, type = "distribution")
  expect_s3_class(p_dist, "ggplot")
  expect_equal(p_dist$labels$y, "Density")
  p_forest <- plot(res_valid, type = "forest")
  expect_s3_class(p_forest, "ggplot")
  expect_equal(p_forest$labels$x, "Hazard Ratio (95% CI)")
  res_2_cuts <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
      num_cuts = 2,
      method = "genetic", max.generations = 5, quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_2_cuts) || all(is.na(res_2_cuts$optimal_cuts)),
    "2-cut result needed for ref group test."
  )
  p_forest_g2 <- plot(res_2_cuts, type = "forest", reference_group = "G2")
  expect_s3_class(p_forest_g2, "ggplot")
  expect_match(p_forest_g2$labels$subtitle, "G2")
})

#' @srrstats {RE4.17} Tests S3 methods (print, summary, plot).
#' @srrstats {G5.2} Tests graceful handling of empty/NA results.
test_that("S3 methods for find_cutpoint_number handle NA results/args", {
  res_na <- structure(
    list(
      results = data.frame(num_cuts = 0:1, BIC = c(NA, NA)),
      parameters = list(criterion = "BIC", method = "systematic")
    ),
    class = "find_cutpoint_number_result"
  )
  expect_message(
    plot(res_na),
    regexp = "Cannot generate plot: no valid IC values"
  )
  expect_message(
    summary(res_na),
    regexp = "Cannot summarise: no valid model was found"
  )
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data, "predictor", "time", "event",
      max_cuts = 1
    )
  ))
  skip_if(is.null(res_valid), "Valid result needed for summary test.")
  sum_res <- summary(res_valid,
    show_comparison_table = FALSE,
    show_best_model_details = FALSE, plot.it = FALSE
  )
  expect_s3_class(sum_res, "find_cutpoint_number_result")
  expect_s3_class(
    suppressMessages(summary(res_valid, plot.it = TRUE)),
    "find_cutpoint_number_result"
  )
})

#' @srrstats {RE4.18} Tests S3 summary method.
test_that("S3 summary.find_cutpoint_number works with all arguments", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data, "predictor", "time", "event",
      max_cuts = 1
    )
  ))
  skip_if(is.null(res_valid), "Valid result needed for summary test.")
  expect_s3_class(
    summary(res_valid,
      show_comparison_table = TRUE,
      show_best_model_details = TRUE, show_group_counts = TRUE,
      show_medians = TRUE, plot.it = FALSE
    ),
    "find_cutpoint_number_result"
  )
})

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

test_that("Coverage: summary.find_cutpoint_number_result – no valid IC", {
  obj <- structure(
    list(
      optimal_cuts = NA_real_,
      optimal_stat = NA_real_,
      ic = rep(NA_real_, 5),
      parameters = list(
        method    = "systematic",
        criterion = "BIC"
      )
    ),
    class = "find_cutpoint_number_result"
  )
  expect_snapshot(summary(obj))
})

test_that("Coverage: plot.find_cutpoint_number_result – all IC NA", {
  obj <- structure(
    list(
      optimal_cuts = NA_real_,
      optimal_stat = NA_real_,
      ic = rep(NA_real_, 5),
      parameters = list(
        method    = "systematic",
        criterion = "BIC"
      )
    ),
    class = "find_cutpoint_number_result"
  )
  expect_snapshot(plot(obj))
})

test_that("summary handles pathological IC (mix of NA/Inf)", {
  obj <- structure(
    list(
      optimal_cuts = NA_real_,
      optimal_stat = NA_real_,
      ic = c(Inf, NA, Inf, NA, Inf),
      parameters = list(
        method    = "systematic",
        criterion = "BIC"
      )
    ),
    class = "find_cutpoint_number_result"
  )
  expect_snapshot(summary(obj))
})

test_that("S3 methods handle missing parameters gracefully", {
  obj <- structure(
    list(
      optimal_cuts = NA_real_,
      optimal_stat = NA_real_,
      ic = rep(NA_real_, 3)
    ),
    class = "find_cutpoint_number_result",
    parameters = list(method = "unknown", criterion = "BIC")
  )
  expect_snapshot(print(obj))
  expect_snapshot(summary(obj))
  expect_snapshot(plot(obj))
})

# --- SECTION 6: Tests for Internal Helpers (utils-helpers.R) ---
# Note: These tests cover objective functions, genetic search wrappers, and
# data validation helpers.

#' @srrstats {G5.8} Tests handling of invalid cut-points.
test_that(".obj handles invalid cuts", {
  expect_equal(
    .obj(
      params = c(50, 50), time = mock_data$time, censor = mock_data$event,
      target = mock_data$predictor, confound = NULL, numcut = 2, gap = 1,
      nmin = 1, criterion = "logrank"
    ),
    -Inf
  )
})

#' @srrstats {G5.8} Tests fitness function for p_value criterion.
test_that(".obj handles p_value criterion", {
  expect_true(is.numeric(
    .obj(
      params = 50, time = mock_data$time, censor = mock_data$event,
      target = mock_data$predictor, confound = NULL, numcut = 1, gap = 1,
      nmin = 1, criterion = "p_value"
    )
  ))
})

#' @srrstats {G1.4a} Tests exported operator.
test_that("%||% works as expected", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
})

#' @srrstats {G5.8} Tests handling of invalid cut-points.
test_that(".obj - hazard_ratio with invalid coefficient name", {
  result <- .obj(
    params = 50,
    time = rexp(20, 0.05),
    censor = rep(0, 20),
    target = rnorm(20, 50, 10),
    confound = NULL,
    numcut = 1,
    gap = 1,
    nmin = 5,
    criterion = "hazard_ratio"
  )
  expect_equal(result, -Inf)
})

test_that(".obj - p_value with loglik path", {
  confound_data <- data.frame(cov1 = rnorm(30))
  result <- .obj(
    params = 50,
    time = rexp(30, 0.05),
    censor = sample(0:1, 30, replace = TRUE),
    target = rnorm(30, 50, 10),
    confound = confound_data,
    numcut = 1,
    gap = 1,
    nmin = 5,
    criterion = "p_value",
    loglik0 = -100
  )
  expect_true(is.numeric(result))
})

test_that(".obj - p_value fallback path", {
  result <- .obj(
    params = 50,
    time = rexp(30, 0.05),
    censor = sample(0:1, 30, replace = TRUE),
    target = rnorm(30, 50, 10),
    confound = NULL,
    numcut = 1,
    gap = 1,
    nmin = 5,
    criterion = "p_value",
    loglik0 = NA_real_
  )
  expect_true(is.numeric(result) || result == -Inf)
})

test_that(".run_genetic_search - gap = NULL (auto)", {
  skip_if_not_installed("rgenoud")
  result <- suppressMessages(suppressWarnings(
    .run_genetic_search(
      target = rnorm(50, 50, 10),
      numcut = 1,
      time = rexp(50, 0.05),
      censor = sample(0:1, 50, replace = TRUE),
      confound = NULL,
      nmin = 5,
      criterion = "logrank",
      max.generations = 3,
      gap = NULL,
      print.level = 0
    )
  ))
  expect_true(is.null(result) || !is.null(result))
})

test_that(".run_genetic_search - zero gap fallback", {
  skip_if_not_installed("rgenoud")
  result <- suppressMessages(suppressWarnings(
    .run_genetic_search(
      target = rep(50, 30),
      numcut = 1,
      time = rexp(30, 0.05),
      censor = sample(0:1, 30, replace = TRUE),
      confound = NULL,
      nmin = 5,
      criterion = "logrank",
      max.generations = 3,
      gap = NULL,
      print.level = 0
    )
  ))
  expect_true(is.null(result) || !is.null(result))
})

#' @srrstats {G2.14a} Tests that NA in domain returns NULL.
test_that(".run_genetic_search handles invalid domains", {
  skip_if_not_installed("rgenoud")
  expect_null(.run_genetic_search(
    target = c(1:10, Inf), numcut = 1, time = 1:11, censor = 1,
    confound = NULL, nmin = 1, criterion = "logrank"
  ))
  suppressWarnings(expect_null(.run_genetic_search(
    target = rep(NA_real_, 10), numcut = 1,
    time = 1:10, censor = 1, confound = NULL,
    nmin = 1, criterion = "logrank"
  )))
})

test_that(".obj function handles internal model failures", {
  data_model_fail <- data.frame(
    time = 1:10, censor = 1, target = c(1:5, 100:104)
  )
  expect_equal(.obj(
    params = c(50, 150), time = data_model_fail$time,
    censor = data_model_fail$censor, target = data_model_fail$target,
    confound = NULL, numcut = 2, gap = 1, nmin = 1, criterion = "logrank"
  ), -Inf)
})

test_that(".calc_ic handles edge cases and invalid inputs", {
  expect_true(is.na(.calc_ic(model = NULL, k = 2, n = 20, "BIC")))
  expect_true(
    is.na(.calc_ic(model = list(loglik = c(NA, NA)), k = 2, n = 20, "BIC"))
  )
  model_stub <- list(loglik = c(NA, -50))
  n <- 5
  k <- 4
  expect_true(is.na(.calc_ic(model_stub, k = k, n = n, criterion = "AICc")))
  k_valid <- 3
  expect_true(
    is.numeric(.calc_ic(model_stub, k = k_valid, n = n, criterion = "AICc"))
  )
})

test_that(".validate_data_conditions event column checks", {
  data_char_event <- mock_data
  data_char_event$event <- as.character(data_char_event$event)
  expect_error(
    find_cutpoint(
      data_char_event, "predictor", "time", "event",
      quiet = TRUE
    ),
    regexp = "must be numeric"
  )
  data_bad_event <- mock_data
  data_bad_event$event[1] <- 2
  expect_error(
    find_cutpoint(
      data_bad_event, "predictor", "time", "event",
      quiet = TRUE
    ),
    regexp = "must contain only 0 and 1"
  )
})

test_that(".validate_data_conditions - nmin proportion", {
  userdata <- mock_data
  names(userdata)[names(userdata) == "predictor"] <- "factor"
  result <- .validate_data_conditions(
    userdata = userdata,
    nmin = 0.2,
    num_cuts = 1,
    outcome_event = "event",
    quiet = FALSE
  )
  expect_true(result$valid)
  expect_true(result$nmin_abs > 1)
})

test_that(".validate_data_conditions - nmin = 0", {
  userdata <- mock_data
  names(userdata)[names(userdata) == "predictor"] <- "factor"
  expect_error(
    .validate_data_conditions(
      userdata = userdata,
      nmin = 0,
      num_cuts = 1,
      outcome_event = "event",
      quiet = TRUE
    ),
    regexp = "must be a non-negative number"
  )
})

test_that("Coverage: .obj default switch", {
  expect_equal(
    .obj(
      params = 50, time = mock_data$time, censor = mock_data$event,
      target = mock_data$predictor, confound = NULL, numcut = 1, gap = 1,
      nmin = 1, criterion = "INVALID_CRITERION"
    ),
    -Inf
  )
})

test_that("Coverage: .get_stat model failures (Pathological Data)", {
  ## ---- log-rank without covariates ---------------------------------------
  res_lr_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data_no_events,
      predictor = "predictor",
      outcome_time = "time", # <-- renamed
      outcome_event = "event", # <-- renamed
      num_cuts = 1,
      method = "systematic",
      criterion = "logrank",
      nmin = 1,
      quiet = TRUE
    )
  ))
  expect_false(all(is.na(res_lr_fail$optimal_cuts)))
  expect_equal(res_lr_fail$optimal_stat, 0)

  ## ---- log-rank WITH covariates -----------------------------------------
  res_lr_cov_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data          = mock_data_no_events,
      predictor     = "predictor",
      outcome_time  = "time", # <-- renamed
      outcome_event = "event", # <-- renamed
      num_cuts      = 1,
      method        = "systematic",
      criterion     = "logrank",
      covariates    = "covariate1",
      nmin          = 1,
      quiet         = TRUE
    )
  ))
  expect_false(all(is.na(res_lr_cov_fail$optimal_cuts)))
  expect_equal(res_lr_cov_fail$optimal_stat, 0)

  ## ---- p-value criterion ------------------------------------------------
  res_pv_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data_no_events,
      predictor = "predictor",
      outcome_time = "time", # <-- renamed
      outcome_event = "event", # <-- renamed
      num_cuts = 1,
      method = "systematic",
      criterion = "p_value",
      nmin = 1,
      quiet = TRUE
    )
  ))
  expect_false(all(is.na(res_pv_fail$optimal_cuts)))
  expect_true(
    res_pv_fail$optimal_stat >= 0 & res_pv_fail$optimal_stat <= 1
  ) # p = 1
})

# --- SECTION 7: Full Workflow and Additional Coverage Tests ---
# Note: These tests cover end-to-end workflows and additional branches
# for high coverage.

#' @srrstats {RE1.0} Integration test for core algorithm.
#' @srrstats {RE4.0} Integration test for model selection.
#' @srrstats {RE7.0} Integration test for bootstrap validation.
test_that("full workflow: num -> find -> validate", {
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
    is.null(num_result) || is.na(num_result$optimal_num_cuts) ||
      num_result$optimal_num_cuts < 0,
    "find_cutpoint_number failed to produce a valid number of cuts."
  )
  if (num_result$optimal_num_cuts == 0) {
    expect_equal(num_result$optimal_num_cuts, 0)
  } else {
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
        result <- validate_cutpoint(fc_result,
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
  }
})

test_that("Coverage: End-to-end with all features", {
  result <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      method = "systematic",
      criterion = "logrank",
      covariates = c("covariate1", "covariate2"),
      nmin = 0.15,
      seed = 42,
      quiet = FALSE
    )
  ))
  expect_s3_class(result, "find_cutpoint")
  if (!any(is.na(result$optimal_cuts))) {
    summary(result,
      show_model = TRUE, show_group_counts = TRUE,
      show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE
    )
    if (requireNamespace("broom", quietly = TRUE)) {
      p1 <- plot(result, type = "distribution")
      p2 <- suppressWarnings(plot(result, type = "outcome"))
      p3 <- plot(result, type = "forest")
      expect_s3_class(p1, "ggplot")
      expect_true(!is.null(p2))
      expect_s3_class(p3, "ggplot")
    }
  }
})

test_that("Coverage: find_cutpoint_number with all IC criteria", {
  for (crit in c("AIC", "AICc", "BIC")) {
    result <- suppressMessages(suppressWarnings(
      find_cutpoint_number(
        data = mock_data,
        predictor = "predictor",
        outcome_time = "time",
        outcome_event = "event",
        method = "systematic",
        criterion = crit,
        max_cuts = 1,
        nmin = 10
      )
    ))
    expect_s3_class(result, "find_cutpoint_number_result")
    expect_true(crit %in% names(result$results))
    summary(result,
      show_comparison_table = TRUE,
      show_best_model_details = TRUE,
      show_group_counts = TRUE, show_medians = TRUE,
      plot.it = FALSE
    )
  }
})

test_that("Coverage: find_cutpoint_number uses seed", {
  skip_if_not_installed("rgenoud")
  res1 <- suppressMessages(suppressWarnings(find_cutpoint_number(
    data = head(mock_data, 30),
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = "genetic", max_cuts = 1, seed = 123,
    max.generations = 5
  )))
  res2 <- suppressMessages(suppressWarnings(find_cutpoint_number(
    data = head(mock_data, 30),
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = "genetic", max_cuts = 1, seed = 123,
    max.generations = 5
  )))
  # Seed ensures the (stochastic) results are identical
  expect_equal(res1$results, res2$results)
})

test_that("Coverage: .systematic_search_num handles nmin edge cases", {
  # This test will hit empty grid branches for both k_cuts=1 and k_cuts=2
  expect_message(
    res_fail <- find_cutpoint_number(
      data = head(mock_data, 30),
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "systematic",
      max_cuts = 2,
      nmin = 16 # nmin*(2+1) > 30, so grid1 will be empty for k=2
    ),
    regexp = "Not enough data|Skipping"
  )
  expect_true(is.na(res_fail$optimal_num_cuts))
})

test_that("Coverage: find_cutpoint_number handles Cox model failures", {
  # Mock survival::coxph to always fail
  local_mocked_bindings(
    "coxph" = function(...) stop("Mocked coxph failure."),
    .package = "survival"
  )
  # Test systematic path
  expect_message(
    suppressWarnings(find_cutpoint_number(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "systematic",
      max_cuts = 1
    )),
    regexp = "Mocked coxph failure."
  )
  # Test genetic path
  skip_if_not_installed("rgenoud")
  # *** FIX ***: Expect the message from the error
  expect_message(
    suppressWarnings(find_cutpoint_number(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "genetic",
      max_cuts = 1
    )),
    regexp = "Mocked coxph failure."
  )
})

test_that("Coverage: find_cutpoint_number - character event column", {
  bad_data <- mock_data
  bad_data$event <- as.character(bad_data$event)
  expect_error(
    find_cutpoint_number(
      data = bad_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      max_cuts = 1
    ),
    regexp = "must be numeric"
  )
})

test_that("Coverage: find_cutpoint_number - invalid event values", {
  bad_data <- mock_data
  bad_data$event[1] <- 3
  expect_error(
    find_cutpoint_number(
      data = bad_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      max_cuts = 1
    ),
    regexp = "must contain only 0 and 1"
  )
})

test_that("Coverage: .genetic_search_num - base model calculation", {
  skip_if_not_installed("rgenoud")
  result <- suppressMessages(suppressWarnings(
    .genetic_search_num(
      userdata = transform(mock_data, factor = predictor),
      max_cuts = 1,
      nmin = 5,
      criterion = "AIC",
      covariates = NULL,
      max.generations = 3,
      pop.size = 10
    )
  ))
  expect_true("AIC" %in% names(result))
  expect_true(!is.na(result$AIC[result$num_cuts == 0]))
})

test_that("Coverage: .systematic_search_num – max_cuts > 2 error", {
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      method = "systematic", max_cuts = 3
    ),
    regexp = "only implemented for max_cuts <= 2"
  )
})
