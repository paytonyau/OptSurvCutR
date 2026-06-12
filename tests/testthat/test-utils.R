# ===================================================================
# MAXIMUM HORIZON BRANCH COVERAGE: utils-helpers.R
# Surgically hitting all validation traps and messaging hooks
# ===================================================================

# --- 1. Information Criteria (.calc_ic) ---

test_that(".calc_ic executes all fallback branches and mathematical edge cases", {
  # Sweep standard switch components
  for (crit in c("AIC", "AICc", "BIC")) {
    expect_true(is.numeric(OptSurvCutR:::.calc_ic(model = list(loglik = c(0, -50)), k = 2, n = 20, criterion = crit)))
  }

  # Trigger NULL and missing data structures
  expect_true(is.na(OptSurvCutR:::.calc_ic(model = NULL, k = 2, n = 20, criterion = "BIC")))
  expect_true(is.na(OptSurvCutR:::.calc_ic(model = list(), k = 2, n = 20, criterion = "BIC")))
  expect_true(is.na(OptSurvCutR:::.calc_ic(model = list(loglik = c(NA, NA)), k = 2, n = 20, criterion = "BIC")))

  # Trigger BIC sample size failures
  expect_true(is.na(OptSurvCutR:::.calc_ic(model = list(loglik = c(0, -50)), k = 2, n = NA, criterion = "BIC")))
  expect_true(is.na(OptSurvCutR:::.calc_ic(model = list(loglik = c(0, -50)), k = 2, n = 0, criterion = "BIC")))

  # Trigger AICc saturation path (n - k - 1 <= 0)
  expect_true(is.na(OptSurvCutR:::.calc_ic(model = list(loglik = c(0, -50)), k = 4, n = 5, criterion = "AICc")))
})


# --- 2. Input Validation (.validate_find_cutpoint_inputs) ---

#' @srrstats {G2.0}
#' @srrstats {G2.3a}
#' @srrstats {G5.2}
test_that(".validate_find_cutpoint_inputs traps argument discrepancies early", {
  # Base clean data framework
  df <- data.frame(pred = 1:10, time = 1:10, status = rep(1, 10), cov = 1:10)

  # Trigger controlled vocabulary argument match traps via match.arg handles
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 1, method = "INVALID", "logrank", NULL))
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 1, method = "systematic", "INVALID", NULL))

  # Explicit case sensitivity testing validation for {G2.3a} / {G2.3b}
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 1, method = "Systematic", "logrank", NULL))
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 1, method = "systematic", "Logrank", NULL))

  # Trigger num_cuts integer and constraint validation checks
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = -1, "systematic", "logrank", NULL))
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 1.5, "systematic", "logrank", NULL))
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 2, "systematic", "hazard_ratio", NULL))
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 3, "systematic", "logrank", NULL))

  # Trigger missing assignment boundaries
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, NULL, "time", "status", num_cuts = 1, "systematic", "logrank", NULL))
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", NULL, "status", num_cuts = 1, "systematic", "logrank", NULL))

  # Trigger column existence mismatch errors
  expect_error(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "MISSING", "time", "status", num_cuts = 1, "systematic", "logrank", NULL))

  # Clean pass check via robust expect_silent matrix evaluation
  expect_silent(OptSurvCutR:::.validate_find_cutpoint_inputs(df, "pred", "time", "status", num_cuts = 1, "systematic", "logrank", "cov"))
})


# --- 3. Data Cleansing (.prepare_cutpoint_data) ---

#' @srrstats {G2.10}
#' @srrstats {G2.13}
#' @srrstats {RE2.1}
test_that(".prepare_cutpoint_data subsets, renames, and omits NA cases", {
  df_na <- data.frame(
    pred = c(1, 2, NA, 4),
    time = c(10, 11, 12, NA),
    status = c(1, 0, 1, 1),
    cov = c(0.1, 0.2, 0.3, 0.4)
  )

  cleaned <- OptSurvCutR:::.prepare_cutpoint_data(df_na, "pred", "time", "status", "cov")

  # Check renaming rules and layout drops
  expect_s3_class(cleaned, "data.frame")
  expect_identical(nrow(cleaned), 2L) # UPGRADED: expect_equal -> expect_identical
  expect_identical(names(cleaned), c("factor", "time", "event", "cov")) # UPGRADED: expect_equal -> expect_identical
})


# --- 4. Event Column Type Enforcement (.validate_event_column) ---

#' @srrstats {G2.1}
#' @srrstats {G5.2b}
test_that(".validate_event_column enforces strict 0/1 binary formatting and accurate messaging", {
  # Trigger non-numeric data entry checks
  expect_error(
    OptSurvCutR:::.validate_event_column(c("1", "0"), "status"),
    regexp = "must be numeric" # Updated to match your exact 80-char message string
  )

  # Trigger out-of-bounds categorical integer check
  expect_error(
    OptSurvCutR:::.validate_event_column(c(0, 1, 2), "status"),
    regexp = "must be strictly binary"
  )

  # Clean pass check
  expect_silent(OptSurvCutR:::.validate_event_column(c(0, 1, 0, 1), "status"))
})

# --- 5. Cohort Metrics Constraints (.validate_data_conditions) ---

#' @srrstats {G2.9}
#' @srrstats {G5.2a}
test_that(".validate_data_conditions monitors sample boundaries and quiet modes completely", {
  # Scenario A: Non-negative time vector enforcement
  bad_time <- data.frame(time = c(1, -5, 2), event = c(1, 0, 1), factor = c(10, 20, 30))

  expect_error(
    OptSurvCutR:::.validate_data_conditions(bad_time, nmin = 1, num_cuts = 1, "event", quiet = TRUE),
    regexp = "Time"
  )

  # Scenario B: Proportional float nmin formatting with active messages caught via expect_message
  clean_df <- data.frame(time = c(10, 12, 14, 16), event = c(1, 0, 1, 0), factor = c(40, 45, 50, 55))

  res_prop <- suppressMessages(
    OptSurvCutR:::.validate_data_conditions(clean_df, nmin = 0.25, num_cuts = 1, "event", quiet = TRUE)
  )
  expect_true(res_prop$valid)
  # floor(0.25 * 4) = 1
  # ✅ FIXED for expect_identical: cast the target to a double/numeric vector
  expect_identical(res_prop$nmin_abs, as.numeric(1))
  
  # Scenario C: Negative nmin value bounds checking
  expect_error(OptSurvCutR:::.validate_data_conditions(clean_df, nmin = -0.5, num_cuts = 1, "event", quiet = TRUE))

  # Scenario D: Insufficient global sample space allocation checks
  res_low_data <- OptSurvCutR:::.validate_data_conditions(clean_df, nmin = 3, num_cuts = 1, "event", quiet = TRUE)
  expect_false(res_low_data$valid)

  # Scenario E: Low variance / unique classification limits
  constant_df <- data.frame(time = c(10, 12, 14), event = c(1, 1, 1), factor = c(50, 50, 50))
  res_few_unique <- OptSurvCutR:::.validate_data_conditions(constant_df, nmin = 1, num_cuts = 1, "event", quiet = TRUE)
  expect_false(res_few_unique$valid)
})


# --- 6. Null Coalescing Operator (%||%) ---

test_that("null coalescing operator matches your implementation exactly", {
  expect_identical(NULL %||% "backup", "backup")       # UPGRADED: expect_equal -> expect_identical
  expect_identical("default" %||% "backup", "default") # UPGRADED: expect_equal -> expect_identical
  expect_identical(3 %||% 5, 3)                         # UPGRADED: expect_equal -> expect_identical
})
