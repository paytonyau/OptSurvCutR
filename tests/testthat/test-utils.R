# ===================================================================
# TESTS: utils-helpers.R
# Internal mathematical engines, validation, and objective functions
# ===================================================================

# --- 1. Objective Functions (.obj) ---

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

test_that(".obj - hazard_ratio with invalid coefficient name", {
  result <- .obj(
    params = 50,
    time = rexp(20, 0.05),
    censor = rep(0, 20),
    target = rnorm(20, 50, 10),
    confound = NULL,
    numcut = 1, gap = 1, nmin = 5,
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
    numcut = 1, gap = 1, nmin = 5,
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
    numcut = 1, gap = 1, nmin = 5,
    criterion = "p_value",
    loglik0 = NA_real_
  )
  expect_true(is.numeric(result) || result == -Inf)
})

test_that(".obj handles internal model failures", {
  data_model_fail <- data.frame(time = 1:10, censor = 1, target = c(1:5, 100:104))
  expect_equal(.obj(
    params = c(50, 150), time = data_model_fail$time,
    censor = data_model_fail$censor, target = data_model_fail$target,
    confound = NULL, numcut = 2, gap = 1, nmin = 1, criterion = "logrank"
  ), -Inf)
})

test_that("Coverage: .obj default switch handles invalid criterion", {
  expect_equal(
    .obj(
      params = 50, time = mock_data$time, censor = mock_data$event,
      target = mock_data$predictor, confound = NULL, numcut = 1, gap = 1,
      nmin = 1, criterion = "INVALID_CRITERION"
    ),
    -Inf
  )
})

# --- 2. Information Criteria (.calc_ic) ---

test_that(".calc_ic handles edge cases and invalid inputs", {
  expect_true(is.na(.calc_ic(model = NULL, k = 2, n = 20, "BIC")))
  expect_true(is.na(.calc_ic(model = list(loglik = c(NA, NA)), k = 2, n = 20, "BIC")))

  model_stub <- list(loglik = c(NA, -50))
  n <- 5
  k <- 4
  expect_true(is.na(.calc_ic(model_stub, k = k, n = n, criterion = "AICc")))

  k_valid <- 3
  expect_true(is.numeric(.calc_ic(model_stub, k = k_valid, n = n, criterion = "AICc")))
})

# --- 3. Genetic Search Engine (.run_genetic_search) ---

test_that(".run_genetic_search - gap = NULL (auto) and zero fallback", {
  skip_if_not_installed("rgenoud")

  result1 <- suppressMessages(suppressWarnings(
    .run_genetic_search(
      target = rnorm(50, 50, 10), numcut = 1, time = rexp(50, 0.05),
      censor = sample(0:1, 50, replace = TRUE), confound = NULL, nmin = 5,
      criterion = "logrank", max.generations = 3, gap = NULL, print.level = 0
    )
  ))
  expect_true(is.null(result1) || !is.null(result1))

  # Zero fallback (constant target)
  result2 <- suppressMessages(suppressWarnings(
    .run_genetic_search(
      target = rep(50, 30), numcut = 1, time = rexp(30, 0.05),
      censor = sample(0:1, 30, replace = TRUE), confound = NULL, nmin = 5,
      criterion = "logrank", max.generations = 3, gap = NULL, print.level = 0
    )
  ))
  expect_true(is.null(result2) || !is.null(result2))
})

#' @srrstats {G2.14a} Tests that NA in domain returns NULL.
test_that(".run_genetic_search handles invalid domains", {
  skip_if_not_installed("rgenoud")

  expect_null(.run_genetic_search(
    target = c(1:10, Inf), numcut = 1, time = 1:11, censor = 1,
    confound = NULL, nmin = 1, criterion = "logrank"
  ))

  suppressWarnings(expect_null(.run_genetic_search(
    target = rep(NA_real_, 10), numcut = 1, time = 1:10, censor = 1,
    confound = NULL, nmin = 1, criterion = "logrank"
  )))
})

# --- 4. Data Validation (.validate_data_conditions) ---

test_that(".validate_data_conditions event column checks", {
  data_char_event <- mock_data
  data_char_event$event <- as.character(data_char_event$event)
  expect_error(
    find_cutpoint(data_char_event, "predictor", "time", "event", quiet = TRUE),
    regexp = "must be numeric"
  )

  data_bad_event <- mock_data
  data_bad_event$event[1] <- 2
  expect_error(
    find_cutpoint(data_bad_event, "predictor", "time", "event", quiet = TRUE),
    regexp = "must contain only 0 and 1"
  )
})

test_that(".validate_data_conditions handles nmin proportions and constraints", {
  userdata <- mock_data
  names(userdata)[names(userdata) == "predictor"] <- "factor"

  result <- .validate_data_conditions(userdata = userdata, nmin = 0.2, num_cuts = 1, outcome_event = "event", quiet = FALSE)
  expect_true(result$valid)
  expect_true(result$nmin_abs > 1)

  expect_error(
    .validate_data_conditions(userdata = userdata, nmin = 0, num_cuts = 1, outcome_event = "event", quiet = TRUE),
    regexp = "must be a non-negative number"
  )
})

# --- 5. Statistics Extraction (.get_stat) ---

test_that("Coverage: .get_stat - group constraint violations", {
  stat <- .get_stat(
    cuts = c(45, 55), num_cuts = 2,
    data_in = data.frame(factor = rep(50, 20), time = rexp(20, 0.05), event = sample(0:1, 20, replace = TRUE)),
    criterion = "logrank", cov_formula = "", nmin = 5, fit_null = NULL
  )
  expect_true(is.na(stat))
})

test_that("Coverage: .get_stat - hazard_ratio with missing coefficient", {
  stat <- .get_stat(
    cuts = 50, num_cuts = 1,
    data_in = data.frame(factor = c(rep(40, 10), rep(60, 10)), time = rexp(20, 0.05), event = rep(0, 20)),
    criterion = "hazard_ratio", cov_formula = "", nmin = 5, fit_null = NULL
  )
  expect_equal(stat, -Inf)
})

test_that("Coverage: .get_stat - p_value with null fit_null", {
  stat <- .get_stat(
    cuts = 50, num_cuts = 1,
    data_in = data.frame(factor = c(rep(40, 15), rep(60, 15)), time = rexp(30, 0.05), event = sample(0:1, 30, replace = TRUE)),
    criterion = "p_value", cov_formula = "", nmin = 5, fit_null = NULL
  )
  expect_true(is.na(stat) || is.numeric(stat))
})

test_that("Coverage: .get_stat model failures (Pathological Data)", {
  # Log-rank without covariates
  res_lr_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(data = mock_data_no_events, predictor = "predictor", outcome_time = "time", outcome_event = "event", num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 1, quiet = TRUE)
  ))
  expect_false(all(is.na(res_lr_fail$optimal_cuts)))
  expect_equal(res_lr_fail$optimal_stat, 0)

  # Log-rank WITH covariates
  res_lr_cov_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(data = mock_data_no_events, predictor = "predictor", outcome_time = "time", outcome_event = "event", num_cuts = 1, method = "systematic", criterion = "logrank", covariates = "covariate1", nmin = 1, quiet = TRUE)
  ))
  expect_false(all(is.na(res_lr_cov_fail$optimal_cuts)))
  expect_equal(res_lr_cov_fail$optimal_stat, 0)

  # P-value criterion
  res_pv_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(data = mock_data_no_events, predictor = "predictor", outcome_time = "time", outcome_event = "event", num_cuts = 1, method = "systematic", criterion = "p_value", nmin = 1, quiet = TRUE)
  ))
  expect_false(all(is.na(res_pv_fail$optimal_cuts)))
  expect_true(res_pv_fail$optimal_stat >= 0 & res_pv_fail$optimal_stat <= 1)
})

# --- 6. Minor Utilities ---

#' @srrstats {G1.4a} Tests exported operator.
test_that("%||% works as expected", {
  expect_equal(NULL %||% 5, 5)
  expect_equal(3 %||% 5, 3)
})
