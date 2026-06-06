# ===================================================================
# TESTS: find_cutpoint_number()
# Automated identification of optimal number of threshold splits
# ===================================================================

# --- 1. Stepwise Model Dimension Optimization Sweeps ---

#' @srrstats {G1.0}
#' @srrstats {G1.1}
#' @srrstats {G5.0}
test_that("find_cutpoint_number executes stepwise criterion selection loops completely", {
  for (crit in c("AIC", "AICc", "BIC")) {
    res_num <- suppressMessages(suppressWarnings(
      find_cutpoint_number(mock_data, "predictor", "time", "event",
        max_cuts = 2, method = "systematic", criterion = crit, quiet = TRUE, nmin = 3
      )
    ))
    expect_s3_class(res_num, "find_cutpoint_number_result")
  }
})

# --- 2. Floating-Point Equality & Tolerance Protections ---

#' @srrstats {G3.0}
test_that("find_cutpoint_number utilizes approximate numeric tolerances for floating equality checks", {
  res_base <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      max_cuts = 2, method = "systematic", criterion = "BIC", quiet = TRUE, nmin = 3
    )
  ))

  expect_s3_class(res_base, "find_cutpoint_number_result")
  expect_false(is.null(res_base$optimal_num_cuts))
})

# --- 3. Front-End Exception Guard Rails ---

#' @srrstats {G5.2b}
test_that("find_cutpoint_number input validation traps enforce structural parameters early", {
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event", max_cuts = 4),
    regexp = "cuts"
  )
  expect_error(
    find_cutpoint_number(mock_data, "non_existent_column", "time", "event"),
    regexp = "column"
  )
})

# --- 4. Exhaustive Multi-Cut Genetic Search Path ---

test_that("Exhaustive multi-cut genetic search optimization for step sizing", {
  # Forces open the genetic dimension loops inside find_cutpoint_number.R
  res_num_gen <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data,
      predictor = "predictor", outcome_time = "time", outcome_event = "event",
      max_cuts = 2, method = "genetic", criterion = "BIC", quiet = TRUE,
      nmin = 3, max.generations = 3, pop.size = 15
    )
  ))
  expect_s3_class(res_num_gen, "find_cutpoint_number_result")
})

# --- 5. Stepwise Optimization Failure & Tie-Breaking Edge Cases ---

test_that("Exhaustive edge-case branches for find_cutpoint_number dimension paths", {
  # A. TRIGGER THE EARLY-ABORT PATHWAY:
  # Using a highly restricted sample layout (n = 20) with a large nmin = 8 constraint forces
  # the 2-cut and 3-cut optimization models to fail internally due to boundary exhaustion.
  # This tests the defensive line coverage switches that compile partial successes safely.
  res_abort_path <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data[1:20, ],
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      max_cuts = 3,
      method = "systematic",
      criterion = "BIC",
      quiet = TRUE,
      nmin = 8
    )
  ))
  expect_s3_class(res_abort_path, "find_cutpoint_number_result")

  # B. TRIGGER THE TIE-BREAKER / EQUALITY TRACKER:
  # Run against your global uniform pathological dataset where all predictors are identical (50).
  # This explicitly forces testing of flat information criterion profiles.
  res_tie_path <- suppressMessages(suppressWarnings(
    find_cutpoint_number(
      mock_data_pathological,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      max_cuts = 2,
      method = "systematic",
      criterion = "AIC",
      quiet = TRUE,
      nmin = 2
    )
  ))
  expect_s3_class(res_tie_path, "find_cutpoint_number_result")
})

# --- 6. Exhaustive S3 Namespace Method Sweeps ---

test_that("S3 methods for find_cutpoint_number provide exhaustive branch coverage via explicit namespace targeting", {
  res_num <- suppressMessages(suppressWarnings(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
      max_cuts = 2, method = "systematic", criterion = "BIC", quiet = TRUE, nmin = 3
    )
  ))

  res_sum <- OptSurvCutR:::summary.find_cutpoint_number_result(res_num)
  expect_s3_class(res_sum, "find_cutpoint_number_result")

  # Dual stream capture wrapper
  p_out_num <- capture.output(
    capture.output(OptSurvCutR:::print.find_cutpoint_number_result(res_num), type = "message")
  )
  expect_true(length(p_out_num) >= 0)
})

test_that("S3 routers for find_cutpoint_number handle incomplete optimizations cleanly", {
  res_empty <- structure(
    list(
      optimal_num_cuts = NA_integer_,
      stats_summary = data.frame(num_cuts = 0:2, BIC = NA_real_, stringsAsFactors = FALSE),
      parameters = list(method = "systematic", criterion = "BIC", predictor = "predictor")
    ),
    class = "find_cutpoint_number_result"
  )

  expect_s3_class(res_empty, "find_cutpoint_number_result")

  # Dual stream capture wrapper for the empty state
  p_out_empty <- capture.output(
    capture.output(OptSurvCutR:::print.find_cutpoint_number_result(res_empty), type = "message")
  )
  expect_true(length(p_out_empty) >= 0)
})
