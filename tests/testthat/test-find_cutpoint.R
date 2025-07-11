# --- Setup: Create a simple dataset for basic tests ---
test_data_binary <- data.frame(
  outcome = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
  predictor = c(1.1, 1.2, 0.9, 1.0, 2.1, 2.0, 1.9, 2.3, 0.8, 2.2)
)

# --- Test Group 1: Input Validation ---
test_that("find_cutpoint stops with incorrect inputs", {
  expect_error(
    find_cutpoint(test_data_binary, predictor = "wrong_name", outcome_binary = "outcome")
  )
})

# --- Test Group 2: Systematic Method ---
test_that("systematic method works for binary outcomes", {
  res_bin <- find_cutpoint(
    data = test_data_binary,
    predictor = "predictor",
    outcome_binary = "outcome",
    nmin = 2
  )
  expect_s3_class(res_bin, "find_cutpoint_systematic")
  expect_true("best_by_vote" %in% names(res_bin))
})


# --- Test Group 3: Genetic Algorithm Method (using a real dataset) ---
test_that("genetic method works for survival data", {
  # The genetic algorithm requires a larger, more realistic dataset to run reliably.
  # We use the built-in 'lung' dataset for this test.
  lung_data <- survival::lung
  lung_data <- na.omit(lung_data)
  lung_data$status <- ifelse(lung_data$status == 2, 1, 0)
  
  set.seed(123)
  # The test is that this complex function runs to completion without error.
  expect_no_error(
    find_cutpoint(
      data = lung_data,
      predictor = "wt.loss",
      outcome_time = "time",
      outcome_event = "status",
      method = "genetic",
      num_cuts = 2,
      nmin = 15,
      numgen = 5 # Use low numgen for a fast test
    )
  )
})