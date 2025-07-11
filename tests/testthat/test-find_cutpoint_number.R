# --- Setup: Create a reliable dataset for testing ---
test_data_surv <- data.frame(
  time = c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
  status = c(1, 0, 1, 0, 1, 0, 1, 0, 1, 0),
  predictor = c(1.1, 1.2, 0.9, 1.0, 2.1, 2.0, 1.9, 2.3, 0.8, 2.2)
)

# --- Test Group 1: Systematic Method with BIC and AIC ---
test_that("systematic method works with BIC and AIC", {
  # Test with BIC (default)
  res_bic <- find_cutpoint_number(
    data = test_data_surv, predictor = "predictor",
    outcome_time = "time", outcome_event = "status",
    method = "systematic", criterion = "BIC", max_cuts = 1, nmin = 2
  )
  expect_s3_class(res_bic, "find_cutpoint_number_result")
  expect_true("BIC" %in% names(res_bic$results))

  # Test with AIC
  res_aic <- find_cutpoint_number(
    data = test_data_surv, predictor = "predictor",
    outcome_time = "time", outcome_event = "status",
    method = "systematic", criterion = "AIC", max_cuts = 1, nmin = 2
  )
  expect_s3_class(res_aic, "find_cutpoint_number_result")
  expect_true("Akaike_Weight" %in% names(res_aic$results))
})


# --- Test Group 2: Genetic Method ---
test_that("genetic method runs without error", {
  set.seed(42)
  # The main test is that it completes without failing
  expect_no_error(
    find_cutpoint_number(
      data = test_data_surv, predictor = "predictor",
      outcome_time = "time", outcome_event = "status",
      method = "genetic", max_cuts = 1, nmin = 2, numgen = 5
    )
  )
})
