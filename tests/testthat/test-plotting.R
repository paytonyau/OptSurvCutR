# --- Setup: Create a reliable dataset for testing ---
test_data_binary <- data.frame(
  outcome = c(0, 0, 0, 0, 0, 1, 1, 1, 1, 1),
  predictor = c(1.1, 1.2, 0.9, 1.0, 2.1, 2.0, 1.9, 2.3, 0.8, 2.2)
)

# --- Test Group 1: Plotting functions run without error ---
test_that("plotting functions execute correctly", {
  
  # Get a result object from a binary analysis
  cut_res <- find_cutpoint(
    data = test_data_binary,
    predictor = "predictor",
    outcome_binary = "outcome",
    method = "systematic",
    nmin = 2
  )
  
  # Test that each plot can be generated without throwing an error
  expect_no_error(plot_effect_size(cut_res))
  expect_no_error(plot_waterfall(cut_res))
})