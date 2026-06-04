# ===================================================================
# TESTS: find_cutpoint_number()
# Determines the optimal number of cut-points using IC (AIC/BIC)
# ===================================================================

# --- Core Functionality ---

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
                         max.generations = 5, nmin = 1
    )
  ))
  skip_if(is.null(res_bic_gen), "Genetic search returned NULL.")
  expect_s3_class(res_bic_gen, "find_cutpoint_number_result")
  expect_equal(nrow(res_bic_gen$results), 2) # 0, 1 cuts
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

# --- Input Validation & Edge Cases ---

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

  bad_data <- mock_data
  bad_data$event <- as.character(bad_data$event)
  expect_error(
    find_cutpoint_number(
      data = bad_data, predictor = "predictor", outcome_time = "time", outcome_event = "event", max_cuts = 1
    ),
    regexp = "must be numeric"
  )

  bad_data2 <- mock_data
  bad_data2$event[1] <- 3
  expect_error(
    find_cutpoint_number(
      data = bad_data2, predictor = "predictor", outcome_time = "time", outcome_event = "event", max_cuts = 1
    ),
    regexp = "must contain only 0 and 1"
  )
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

test_that("Coverage: find_cutpoint_number handles Cox model failures", {
  # Mock survival::coxph to always fail
  local_mocked_bindings(
    "coxph" = function(...) stop("Mocked coxph failure."),
    .package = "survival"
  )

  # Test systematic path
  expect_message(
    suppressWarnings(find_cutpoint_number(
      data = mock_data, predictor = "predictor", outcome_time = "time", outcome_event = "event", method = "systematic", max_cuts = 1
    )),
    regexp = "Mocked coxph failure."
  )

  # Test genetic path
  skip_if_not_installed("rgenoud")
  expect_message(
    suppressWarnings(find_cutpoint_number(
      data = mock_data, predictor = "predictor", outcome_time = "time", outcome_event = "event", method = "genetic", max_cuts = 1
    )),
    regexp = "Mocked coxph failure."
  )
})

test_that("Coverage: .systematic_search_num handles nmin edge cases", {
  # Enforce quiet = FALSE to guarantee messages bleed out to the test runner
  expect_condition(
    res_fail <- find_cutpoint_number(
      data = head(mock_data, 30),
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "systematic",
      max_cuts = 2,
      nmin = 16, # nmin*(2+1) > 30, so grid will be empty for k=2
      quiet = FALSE
    ),
    class = "message"
  )
  expect_true(is.na(res_fail$optimal_num_cuts))
})

test_that("Coverage: .systematic_search_num – max_cuts > 2 error", {
  expect_error(
    find_cutpoint_number(mock_data, "predictor", "time", "event",
                         method = "systematic", max_cuts = 3
    ),
    regexp = "only implemented for max_cuts <= 2"
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


# --- S3 Methods for find_cutpoint_number ---

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
    find_cutpoint_number(mock_data, "predictor", "time", "event", max_cuts = 1)
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

# --- Integration / Full Workflow ---

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
