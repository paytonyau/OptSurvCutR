# ===================================================================
# WORKFLOW INTEGRATION TESTS
#
# The three core functions are individually well covered, but both defects
# identified during v0.11.0 development were handoff failures between them
# rather than faults within any one function:
#
#   1. `userdata` was dropped between `validate_cutpoint()` and its summary
#      method, so the stability denominator silently fell back to a wrong
#      quantity.
#   2. `num_cuts` was hard-coded in downstream analysis scripts rather than
#      carried from `find_cutpoint_number()`, allowing the reported model
#      complexity to diverge from the selected complexity.
#
# Unit tests cannot detect this class of failure. These tests run the complete
# workflow and assert that the objects connect.
#
# Note: the shared `mock_data_3groups` fixture generates `time` and `event`
# independently of `predictor`, so no threshold exists to recover and model
# selection correctly returns k = 0. These tests therefore generate data with a
# genuine survival signal locally.
# ===================================================================

#' @srrstats {G5.4} Verifies the workflow's internal consistency end to end.
#' @srrstats {G5.5} Fixed seeds throughout for deterministic comparison.
#' @srrstats {RE4.0} Confirms each stage returns the documented class.
#' @noRd
NULL

# Two-threshold data: hazard rises in steps at 40 and 70.
make_workflow_data <- function(n = 300, seed = 1, with_covariate = TRUE) {
  set.seed(seed)
  predictor <- stats::runif(n, 0, 100)
  step <- ifelse(predictor > 70, 4, ifelse(predictor > 40, 2, 1))
  d <- data.frame(
    predictor = predictor,
    time = stats::rexp(n, rate = 0.02 * step),
    event = 1L
  )
  if (with_covariate) d$covariate1 <- stats::rnorm(n, 5, 1)
  d
}

test_that("the three-step workflow preserves state across functions", {
  # skip_on_cran(): fits a 2-cut systematic model (99 candidate positions)
  # inside a 30-replicate bootstrap - observed at several minutes locally.
  # Full coverage of this path is retained via CI and local runs, which is
  # what matters for a handoff-correctness test; CRAN's check-time budget
  # is better spent on the cheaper tests in this file.
  skip_on_cran()

  d <- make_workflow_data(n = 300, seed = 1)

  num <- find_cutpoint_number(
    data = d,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = "systematic",
    criterion = "BIC",
    max_cuts = 2,
    nmin = 0.20,
    seed = 1
  )
  expect_s3_class(num, "find_cutpoint_number_result")

  k <- num$optimal_num_cuts
  skip_if(is.null(k) || is.na(k) || k < 1,
          "No cut-point selected; nothing to hand off")

  cut <- find_cutpoint(
    data = d,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = if (k <= 2) "systematic" else "genetic",
    criterion = "logrank",
    num_cuts = k, # carried forward, not hard-coded
    nmin = 0.20,
    n_perm = 0,
    seed = 1
  )
  expect_s3_class(cut, "find_cutpoint")

  # The number of thresholds located must equal the number selected
  expect_equal(length(cut$optimal_cuts), k)

  val <- validate_cutpoint(cut, num_replicates = 30, n_cores = 1, seed = 1)
  expect_s3_class(val, "validate_cutpoint_result")

  # One confidence interval per located threshold
  expect_equal(nrow(val$confidence_intervals), length(cut$optimal_cuts))

  # The predictor must be the same data at every stage
  expect_identical(val$userdata$factor, cut$userdata$factor)

  # A tier must be assignable from a complete, valid workflow
  s <- summary(val)$stability
  expect_false(is.na(s$tier))
  expect_false(is.na(s$percent))
})

test_that("confidence intervals are ordered consistently with the thresholds", {
  # skip_on_cran(): same cost driver as the previous test - a 2-cut
  # systematic search (99 candidates) inside a bootstrap loop.
  skip_on_cran()

  d <- make_workflow_data(n = 300, seed = 2)

  cut <- find_cutpoint(
    data = d,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = "systematic",
    criterion = "logrank",
    num_cuts = 2,
    nmin = 0.20,
    n_perm = 0,
    seed = 2
  )
  skip_if(anyNA(cut$optimal_cuts), "No valid two-cut solution found")

  val <- validate_cutpoint(cut, num_replicates = 40, n_cores = 1, seed = 2)

  # Thresholds are returned in ascending order, so the intervals must be too
  expect_identical(cut$optimal_cuts, sort(cut$optimal_cuts))
  expect_identical(
    val$confidence_intervals$Lower,
    sort(val$confidence_intervals$Lower)
  )

  # Intervals must be non-degenerate
  widths <- val$confidence_intervals$Upper - val$confidence_intervals$Lower
  expect_true(all(widths >= 0))

  # One relative width is reported per threshold
  s <- summary(val)$stability
  expect_equal(length(s$relative_widths), nrow(val$confidence_intervals))
})

test_that("covariate adjustment is carried through the whole workflow", {
  d <- make_workflow_data(n = 300, seed = 3)
  covs <- "covariate1"

  cut <- find_cutpoint(
    data = d,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    covariates = covs,
    method = "systematic",
    criterion = "logrank",
    num_cuts = 1,
    nmin = 0.20,
    n_perm = 0,
    seed = 3
  )
  skip_if(anyNA(cut$optimal_cuts), "No valid single-cut solution found")

  val <- validate_cutpoint(cut, num_replicates = 30, n_cores = 1, seed = 3)

  # The covariates used in the search must be recorded and reused in
  # validation, otherwise the bootstrap estimates a different model
  expect_identical(val$parameters$covariates, covs)
  expect_true(all(covs %in% names(val$userdata)))
})

test_that("num_cuts = 0 is rejected rather than fitted", {
  # `find_cutpoint_number()` may legitimately return k = 0 for a predictor with
  # no relationship to the outcome. Passing that value onward must fail
  # clearly rather than producing an uninterpretable result.
  d <- make_workflow_data(n = 120, seed = 4)

  expect_error(
    find_cutpoint(d, "predictor", "time", "event",
      num_cuts = 0, method = "systematic",
      nmin = 0.20, n_perm = 0, seed = 4
    )
  )
})
