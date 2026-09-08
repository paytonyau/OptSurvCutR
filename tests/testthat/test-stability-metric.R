# ===================================================================
# STABILITY METRIC REGRESSION TESTS
#
# Guards the defect fixed in v0.11.0: `validate_cutpoint()` did not carry
# `userdata` into its returned object, so `summary()` computed the relative
# confidence interval width against the bootstrap distribution of the first
# cut-point instead of the predictor's 10th-90th percentile range. Widths were
# inflated and Stability Tier 1 was unreachable for every input.
#
# The reachability test below is the one that would have caught it: the failure
# mode was an entire classification outcome that no code path could produce,
# which value-comparison tests cannot detect.
# ===================================================================

#' @srrstats {G5.2} Asserts informative failure when the predictor is unavailable.
#' @srrstats {G5.8} Covers the degenerate case of an undefined percentile range.
#' @srrstats {G5.9} The stability denominator depends only on the observed
#'   predictor and is therefore invariant to the bootstrap seed.
#' @noRd
NULL

# Predictor with a genuine threshold, large enough for stable intervals.
# Generated locally rather than in helper-data.R because the reachability test
# requires a wider spread and more events than the shared mock data provides.
make_threshold_data <- function(n = 800, seed = 1, cut = 50, hr = 3) {
  set.seed(seed)
  marker <- stats::runif(n, 0, 100)
  rate <- 0.02 * ifelse(marker > cut, hr, 1)
  data.frame(
    marker = marker,
    time = stats::rexp(n, rate = rate),
    event = 1L
  )
}

fit_single_cut <- function(d, seed) {
  find_cutpoint(d, "marker", "time", "event",
    num_cuts = 1, method = "systematic",
    nmin = 0.25, n_perm = 0, seed = seed
  )
}

# ------------------------------------------------------------------
# 1. The predictor must survive into the validation object
# ------------------------------------------------------------------

test_that("validate_cutpoint() carries userdata into its result", {
  d <- make_threshold_data(n = 300, seed = 11)
  cp <- fit_single_cut(d, 11)
  val <- validate_cutpoint(cp, num_replicates = 30, n_cores = 1, seed = 11)

  expect_false(is.null(val$userdata))
  expect_true("factor" %in% names(val$userdata))
  expect_identical(val$userdata$factor, cp$userdata$factor)
})

# ------------------------------------------------------------------
# 2. The denominator must be the predictor, not the bootstrap spread
# ------------------------------------------------------------------

test_that("relative width is computed against the predictor's P90-P10 range", {
  d <- make_threshold_data(n = 300, seed = 12)
  cp <- fit_single_cut(d, 12)
  val <- validate_cutpoint(cp, num_replicates = 40, n_cores = 1, seed = 12)
  s <- summary(val)$stability

  expected_spread <- as.numeric(diff(stats::quantile(d$marker, c(0.10, 0.90))))
  expect_equal(s$data_spread, expected_spread, tolerance = 1e-8)

  width <- max(val$confidence_intervals$Upper - val$confidence_intervals$Lower)
  expect_equal(s$percent, round(100 * width / expected_spread, 1), tolerance = 0.2)

  # Explicit guard against the regression: the bootstrap distribution of a
  # cut-point must not be used as the denominator.
  boot_spread <- as.numeric(diff(stats::quantile(
    val$bootstrap_distribution[, 1], c(0.10, 0.90), na.rm = TRUE
  )))
  if (!isTRUE(all.equal(boot_spread, expected_spread))) {
    expect_false(isTRUE(all.equal(s$data_spread, boot_spread)))
  }
})

# ------------------------------------------------------------------
# 3. Tier 1 must be attainable
# ------------------------------------------------------------------

test_that("Tier 1 is reachable for a well-separated threshold", {
  # Before v0.11.0 the relative width always exceeded 100%, so no input of any
  # kind could be classified Tier 1.
  d <- make_threshold_data(n = 800, seed = 13, cut = 50, hr = 4)
  cp <- fit_single_cut(d, 13)
  val <- validate_cutpoint(cp, num_replicates = 150, n_cores = 1, seed = 13)
  s <- summary(val)$stability

  expect_identical(s$tier, "Tier 1")
  expect_lt(s$max_relative_width, 0.30)
})

# ------------------------------------------------------------------
# 4. Single cut-point grading
# ------------------------------------------------------------------

test_that("single cut-point models are graded on width alone", {
  # Interval separation is undefined with one boundary, so Tier 3, which
  # requires overlap, must never be assigned.
  d <- make_threshold_data(n = 300, seed = 14, hr = 1.5)
  cp <- fit_single_cut(d, 14)
  val <- validate_cutpoint(cp, num_replicates = 40, n_cores = 1, seed = 14)
  s <- summary(val)$stability

  expect_true(s$tier %in% c("Tier 1", "Tier 2", "Tier 4"))
  expect_false(identical(s$tier, "Tier 3"))
  expect_true(is.na(s$separated))
})

# ------------------------------------------------------------------
# 5. Undefined cases must report no tier rather than a wrong one
# ------------------------------------------------------------------

test_that("a missing predictor yields no tier rather than an incorrect one", {
  d <- make_threshold_data(n = 200, seed = 15)
  cp <- fit_single_cut(d, 15)
  val <- validate_cutpoint(cp, num_replicates = 30, n_cores = 1, seed = 15)

  val$userdata <- NULL # reproduce the pre-0.11.0 state

  s <- suppressMessages(summary(val))$stability
  expect_true(is.na(s$tier))
  expect_true(is.na(s$percent))
})

# ------------------------------------------------------------------
# 6. The returned contract
# ------------------------------------------------------------------

test_that("summary() returns the documented stability fields", {
  # Asserts the public return contract. The v0.11.0 defect arose because an
  # element was dropped from a returned object without any test noticing.
  d <- make_threshold_data(n = 250, seed = 16)
  cp <- fit_single_cut(d, 16)
  val <- validate_cutpoint(cp, num_replicates = 30, n_cores = 1, seed = 16)
  s <- summary(val)

  expect_s3_class(s, "validate_cutpoint_result")
  expect_true(is.list(s$stability))
  expect_named(
    s$stability,
    c(
      "tier", "tier_label", "percent", "max_relative_width",
      "relative_widths", "data_spread", "separated", "worst_cut"
    ),
    ignore.order = TRUE
  )

  # summary() must not alter the rest of the object
  expect_identical(s$original_cuts, val$original_cuts)
  expect_identical(s$confidence_intervals, val$confidence_intervals)
})
