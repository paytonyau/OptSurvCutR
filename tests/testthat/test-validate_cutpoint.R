# --- Test Group 1: Validation works for systematic results ---
test_that("validate_cutpoint works for systematic results", {

  # Create a dummy result object to pass to the validation function
  dummy_systematic_result <- list(
    userdata = data.frame(factor = 1:20, time = 1:20, event = 1:20),
    best_by_vote = 15,
    parameters = list(analysis_type = "survival")
  )
  class(dummy_systematic_result) <- "find_cutpoint_systematic"

  # Create a simple "mock" function that will temporarily replace the real
  # find_cutpoint. This mock function is guaranteed to return a valid result.
  mock_find_cutpoint <- function(...) {
    return(dummy_systematic_result)
  }

  # Use testthat::with_mocked_bindings to run the test.
  testthat::with_mocked_bindings(
    find_cutpoint = mock_find_cutpoint,
    {
      val_res <- validate_cutpoint(
        cutpoint_result = dummy_systematic_result,
        num_replicates = 10,
        use_parallel = FALSE # Force sequential execution for the test
      )

      # The test now checks if the validation function correctly processed the mock results
      expect_s3_class(val_res, "validate_cutpoint_result")
      expect_true("confidence_intervals" %in% names(val_res))
    }
  )
})


# --- Test Group 2: Validation works for genetic results ---
test_that("validate_cutpoint works for genetic results", {

  # Create a dummy result object for the genetic method
  dummy_genetic_result <- list(
    userdata = data.frame(factor = 1:20, time = 1:20, event = 1:20),
    optimal_cuts = c(10, 20), # A result with 2 cut-points
    parameters = list(analysis_type = "survival")
  )
  class(dummy_genetic_result) <- "find_cutpoint_genetic"

  # Create a mock function that always returns the dummy genetic result
  mock_find_cutpoint_gen <- function(...) {
    return(dummy_genetic_result)
  }

  # Run the test using the mock function
  testthat::with_mocked_bindings(
    find_cutpoint = mock_find_cutpoint_gen,
    {
      val_res_gen <- validate_cutpoint(
        cutpoint_result = dummy_genetic_result,
        num_replicates = 10,
        use_parallel = FALSE # Force sequential execution for the test
      )

      expect_s3_class(val_res_gen, "validate_cutpoint_result")
      # Check that it correctly created CIs for two cut-points
      expect_equal(nrow(val_res_gen$confidence_intervals), 2)
    }
  )
})
