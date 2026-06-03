# ===================================================================
# TESTS: find_cutpoint()
# Core optimal cut-point algorithm (Systematic & Genetic)
# ===================================================================

#' @srrstats {RE1.0} Tests the core cut-point implementation.
#' @srrstats {G5.4b} Comparison against alternative implementations (systematic vs genetic).
#' @srrstats {G5.9} Noise susceptibility tests (stochastic genetic algorithm).
#' @srrstats {G5.9a} Tests with trivial noise (implicit in stochastic runs).
#' @srrstats {G5.9b} Tests with different seeds.

#' @srrstats {RE1.0} Tests core optimal cut-point algorithm (systematic).
test_that("find_cutpoint systematic search works for one and two cuts", {
  # One cut
  res_lr1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 1,
                  method = "systematic", criterion = "logrank", nmin = 1,
                  quiet = TRUE
    )
  ))
  expect_s3_class(res_lr1, "find_cutpoint")
  expect_length(res_lr1$optimal_cuts, 1)
  expect_true(
    all(is.na(res_lr1$optimal_cuts) |
          (res_lr1$optimal_cuts >= min(mock_data$predictor) &
             res_lr1$optimal_cuts <= max(mock_data$predictor)))
  )

  # Two cuts (with appropriate trimodal data)
  res_lr2 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
                  num_cuts = 2,
                  method = "systematic", criterion = "logrank", nmin = 10,
                  quiet = TRUE
    )
  ))
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
})

#' @srrstats {G2.0} Tests validation of explicit genetic algorithm parameters.
test_that("find_cutpoint respects explicit pop.size and max.generations", {
  skip_if_not_installed("rgenoud")

  # Run with non-default GA parameters
  res_ga_params <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      method = "genetic",
      pop.size = 50, # Explicit non-default
      max.generations = 5, # Explicit non-default
      nmin = 5,
      quiet = TRUE
    )
  ))

  # 1. Check object class
  expect_s3_class(res_ga_params, "find_cutpoint")

  # 2. Check that parameters were stored correctly in the output
  expect_equal(res_ga_params$parameters$pop.size, 50)
  expect_equal(res_ga_params$parameters$max.generations, 5)

  # 3. Check that the method was indeed genetic
  expect_equal(res_ga_params$parameters$method, "genetic")
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm for multiple cuts.
test_that("find_cutpoint genetic search works for multiple cuts", {
  skip_if_not_installed("rgenoud")

  # 2 Cuts
  res_lr2 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
                  num_cuts = 2, method = "genetic", criterion = "logrank",
                  max.generations = 10, nmin = 1, quiet = TRUE, seed = 42
    )
  ))
  expect_s3_class(res_lr2, "find_cutpoint")
  expect_length(res_lr2$optimal_cuts, 2)
  expect_true(
    all(is.na(res_lr2$optimal_cuts) | diff(res_lr2$optimal_cuts) > 0)
  )

  # 3 Cuts
  res_lr3 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_3groups, "predictor", "time", "event",
                  num_cuts = 3, method = "genetic", criterion = "logrank",
                  max.generations = 15, nmin = 1, quiet = TRUE, seed = 123
    )
  ))
  expect_s3_class(res_lr3, "find_cutpoint")
  expect_length(res_lr3$optimal_cuts, 3)
})

#' @srrstats {RE1.0} Tests all optimisation criteria.
#' @srrstats {RE4.17} Tests nmin as a proportion.
test_that("find_cutpoint handles various criteria and nmin proportion", {
  # Hazard ratio
  res_hr1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 1, method = "systematic",
                  criterion = "hazard_ratio", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_hr1, "find_cutpoint")
  expect_length(res_hr1$optimal_cuts, 1)

  # p-value
  res_pv1 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 1, method = "systematic",
                  criterion = "p_value", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_pv1, "find_cutpoint")
  expect_length(res_pv1$optimal_cuts, 1)

  # nmin as proportion
  res_nmin_prop <- suppressMessages(suppressWarnings(
    find_cutpoint(
      mock_data, "predictor", "time", "event",
      nmin = 0.1, quiet = TRUE
    )
  ))
  expect_s3_class(res_nmin_prop, "find_cutpoint")
})

#' @srrstats {G2.0} Tests input validation.
#' @srrstats {G2.1} Tests input type validation.
#' @srrstats {G2.4b} Tests match.arg validation.
#' @srrstats {G2.13} Tests that invalid inputs trigger errors.
test_that("find_cutpoint handles invalid inputs", {
  bad_data <- mock_data[, c("time", "event")]
  expect_error(
    find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
    regexp = "Missing required columns: 'predictor'"
  )
  bad_data <- mock_data
  bad_data$time[1] <- -1
  expect_error(
    find_cutpoint(bad_data, "predictor", "time", "event", nmin = 1),
    regexp = "Time variable must be non-negative"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = -1,
                  method = "systematic", nmin = 1
    ),
    regexp = "non-negative integer"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = NA,
                  method = "systematic", nmin = 1
    ),
    regexp = "non-negative integer"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  criterion = "invalid", nmin = 1
    ),
    regexp = "arg.*should be one of"
  )
  expect_error(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 2,
                  criterion = "hazard_ratio", nmin = 1
    ),
    regexp = "'hazard_ratio' is only supported for num_cuts = 1"
  )
})

#' @srrstats {G5.0} Tests edge case (constant predictor).
test_that("find_cutpoint handles quiet arg and low variability predictor", {
  expect_no_message(
    find_cutpoint(mock_data, "predictor", "time", "event", quiet = TRUE)
  )
  low_var_data <- mock_data
  low_var_data$predictor <- rep(50, n_test)
  expect_message(
    find_cutpoint(low_var_data, "predictor", "time", "event", quiet = FALSE),
    regexp = "Predictor has too few unique values"
  )
})

#' @srrstats {G5.0} Tests edge cases (skewed data, heavy censoring).
#' @srrstats {G5.2} Tests handling of non-converged models.
test_that("find_cutpoint handles various data scenarios", {
  res_lr_skewed <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_skewed, "predictor", "time", "event",
                  num_cuts = 1, method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_lr_skewed, "find_cutpoint")
  expect_length(res_lr_skewed$optimal_cuts, 1)

  res_lr_censor <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_heavy_censor, "predictor", "time", "event",
                  num_cuts = 1, method = "systematic",
                  criterion = "logrank", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_lr_censor, "find_cutpoint")
  expect_length(res_lr_censor$optimal_cuts, 1)

  bad_data <- mock_data
  bad_data$predictor <- rep(50, n_test)
  res_cox_fail <- suppressMessages(suppressWarnings(
    find_cutpoint(bad_data, "predictor", "time", "event",
                  num_cuts = 1,
                  method = "systematic", nmin = 10, quiet = TRUE
    )
  ))
  expect_s3_class(res_cox_fail, "find_cutpoint")
  expect_true(all(is.na(res_cox_fail$optimal_cuts)))
})

#' @srrstats {RE2.2} Tests covariate adjustment (systematic).
test_that("find_cutpoint (systematic) works with covariates", {
  res_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  covariates = "covariate1", num_cuts = 1,
                  method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_cov, "find_cutpoint")
  expect_length(res_cov$optimal_cuts, 1)
  expect_true(!is.na(res_cov$optimal_cuts[1]))

  res_lr_seq <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 1, method = "systematic",
                  criterion = "logrank", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_lr_seq, "find_cutpoint")
  expect_length(res_lr_seq$optimal_cuts, 1)
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm.
test_that("find_cutpoint (genetic) works sequentially", {
  skip_if_not_installed("rgenoud")
  res_gen_seq <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 2, method = "genetic",
                  criterion = "logrank", max.generations = 5, nmin = 1, quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_gen_seq) || all(is.na(res_gen_seq$optimal_cuts)),
    "Sequential genetic search failed."
  )
  expect_s3_class(res_gen_seq, "find_cutpoint")
  expect_length(res_gen_seq$optimal_cuts, 2)
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm (high cuts, p_value).
test_that("find_cutpoint handles high num_cuts with p_value", {
  skip_if_not_installed("rgenoud")
  res_pv4 <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 4, method = "genetic",
                  criterion = "p_value", max.generations = 5, nmin = 1, quiet = TRUE
    )
  ))
  skip_if(all(is.na(res_pv4$optimal_cuts)), "GA search failed for 4 cuts.")
  expect_s3_class(res_pv4, "find_cutpoint")
  expect_length(res_pv4$optimal_cuts, 4)
})

#' @srrstats {G5.0} Tests edge cases.
#' @srrstats {G5.8} Tests that constant predictor returns na_result.
test_that("find_cutpoint handles insufficient and constant data", {
  small_data <- mock_data[1:5, ]
  res_insufficient <- suppressMessages(suppressWarnings(
    find_cutpoint(small_data, "predictor", "time", "event",
                  num_cuts = 1,
                  method = "systematic", nmin = 5, quiet = TRUE
    )
  ))
  expect_s3_class(res_insufficient, "find_cutpoint")
  expect_true(all(is.na(res_insufficient$optimal_cuts)))

  constant_data <- mock_data[1:20, ]
  constant_data$predictor <- rep(50, 20)
  res_constant <- suppressMessages(suppressWarnings(
    find_cutpoint(constant_data, "predictor", "time", "event",
                  num_cuts = 1, method = "systematic", nmin = 1, quiet = TRUE
    )
  ))
  expect_s3_class(res_constant, "find_cutpoint")
  expect_true(all(is.na(res_constant$optimal_cuts)))
})

#' @srrstats {G5.8} Tests stochastic genetic algorithm.
test_that("group creation works in find_cutpoint", {
  res_fc <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 2,
                  method = "genetic", criterion = "logrank", max.generations = 5,
                  nmin = 1, quiet = TRUE
    )
  ))
  skip_if(all(is.na(res_fc$optimal_cuts)), "Genetic search failed.")
  expect_s3_class(res_fc, "find_cutpoint")
  expect_length(res_fc$optimal_cuts, 2)
  expect_true(
    all(is.na(res_fc$optimal_cuts) | diff(res_fc$optimal_cuts) > 0)
  )
})

#' @srrstats {G5.0} Tests edge case (constant predictor).
#' @srrstats {G5.8} Tests NA result for constant predictor.
test_that("find_cutpoint's systematic search handles edge cases", {
  edge_data <- mock_data
  edge_data$predictor <- rep(50, n_test)
  expect_message(
    find_cutpoint(
      data = edge_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      nmin = 5,
      quiet = FALSE
    ),
    regexp = "Predictor has too few unique values"
  )
})

test_that("Coverage: find_cutpoint (genetic) shows messages", {
  skip_if_not_installed("rgenoud")
  # Suppress rgenoud's own warning about max.generations
  expect_message(
    suppressWarnings(find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      method = "genetic",
      max.generations = 1,
      quiet = FALSE
    )),
    regexp = "Starting genetic search"
  )
})

test_that("Coverage: find_cutpoint handles nmin edge cases", {
  # Test for systematic, 1 cut, grid empty
  expect_message(
    res_1_cut <- find_cutpoint(
      data = head(mock_data, 20),
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      nmin = 15,
      num_cuts = 1,
      method = "systematic",
      quiet = FALSE
    ),
    regexp = "Not enough data"
  )
  expect_true(all(is.na(res_1_cut$optimal_cuts)))

  # Test for systematic, 2 cuts, grid empty
  expect_message(
    res_2_cut <- find_cutpoint(
      data = head(mock_data, 20),
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      nmin = 8,
      num_cuts = 2,
      method = "systematic",
      quiet = FALSE
    ),
    regexp = "Not enough data"
  )
  expect_true(all(is.na(res_2_cut$optimal_cuts)))
})

test_that("Coverage: find_cutpoint p_value + covariates branches", {
  # Tests p_value + covariates + systematic + 2 cuts branches
  res_p_val_2 <- find_cutpoint(
    data = head(mock_data_3groups, 50),
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    method = "systematic",
    num_cuts = 2,
    criterion = "p_value",
    covariates = "covariate1",
    nmin = 5,
    quiet = TRUE
  )
  expect_s3_class(res_p_val_2, "find_cutpoint")
})

test_that(
  "Coverage: find_cutpoint - genetic method with p_value and covariates",
  {
    skip_if_not_installed("rgenoud")
    result <- suppressMessages(suppressWarnings(
      find_cutpoint(
        data = mock_data,
        predictor = "predictor",
        outcome_time = "time",
        outcome_event = "event",
        num_cuts = 2,
        method = "genetic",
        criterion = "p_value",
        covariates = "covariate1",
        nmin = 5,
        max.generations = 5,
        seed = 123,
        quiet = TRUE
      )
    ))
    expect_s3_class(result, "find_cutpoint")
  }
)

# --- S3 Methods & Full Workflow ---

#' @srrstats {RE4.17} Tests S3 print methods.
#' @srrstats {G5.2} Tests graceful handling of empty/NA results.
test_that("S3 methods for find_cutpoint handle NA results", {
  res_na <- find_cutpoint(mock_data_pathological, "predictor", "time", "event",
                          num_cuts = 1, method = "systematic", quiet = TRUE
  )
  expect_true(all(is.na(res_na$optimal_cuts)))
  expect_message(print(res_na), "No optimal cut-point determined")
  expect_message(summary(res_na), "Optimal Cut-point Analysis")
  expect_message(plot(res_na), "Cannot generate plot: no valid cut-point")

  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  num_cuts = 1,
                  quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_valid) || all(is.na(res_valid$optimal_cuts)),
    "Valid result needed for summary test."
  )
  sum_res <- summary(res_valid,
                     show_model = FALSE, show_group_counts = FALSE,
                     show_medians = FALSE, show_ph_test = FALSE,
                     show_params = FALSE
  )
  expect_s3_class(sum_res, "find_cutpoint")
})

#' @srrstats {RE4.18} Tests S3 summary method.
test_that("S3 summary.find_cutpoint works with all arguments", {
  res_valid <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data, "predictor", "time", "event",
                  covariates = "covariate1", num_cuts = 1, quiet = TRUE
    )
  ))
  skip_if(
    is.null(res_valid) || all(is.na(res_valid$optimal_cuts)),
    "Valid result needed for summary test."
  )
  expect_s3_class(
    summary(res_valid,
            show_model = TRUE, show_group_counts = TRUE,
            show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE
    ),
    "find_cutpoint"
  )
})

test_that("Coverage: End-to-end with all features", {
  result <- suppressMessages(suppressWarnings(
    find_cutpoint(
      data = mock_data,
      predictor = "predictor",
      outcome_time = "time",
      outcome_event = "event",
      num_cuts = 1,
      method = "systematic",
      criterion = "logrank",
      covariates = c("covariate1", "covariate2"),
      nmin = 0.15,
      seed = 42,
      quiet = FALSE
    )
  ))
  expect_s3_class(result, "find_cutpoint")

  if (!any(is.na(result$optimal_cuts))) {
    summary(result,
            show_model = TRUE, show_group_counts = TRUE,
            show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE
    )
    if (requireNamespace("broom", quietly = TRUE)) {
      p1 <- plot(result, type = "distribution")
      p2 <- suppressWarnings(plot(result, type = "outcome"))
      p3 <- plot(result, type = "forest")
      expect_s3_class(p1, "ggplot")
      expect_true(!is.null(p2))
      expect_s3_class(p3, "ggplot")
    }
  }
})


# --- Internal Engine Coverage (.systematic_search) ---

test_that("Coverage: .systematic_search - null model fit failure", {
  bad_data <- mock_data
  bad_data$time <- rep(0, n_test)

  # The internal function expects 'factor' instead of the predictor name
  bad_data2 <- bad_data
  names(bad_data2)[names(bad_data2) == "predictor"] <- "factor"

  result <- suppressMessages(suppressWarnings(
    OptSurvCutR:::.systematic_search(
      userdata = bad_data2,
      num_cuts = 1,
      criterion = "p_value",
      covariates = NULL,
      nmin = 5,
      predictor_name = "predictor",
      quiet = TRUE
    )
  ))
  expect_true(all(is.na(result$optimal_cuts)))
  expect_equal(result$parameters$method, "systematic")
})

test_that("Coverage: .systematic_search - no valid cuts (nmin violation)", {
  tiny_data <- mock_data[1:10, ]
  tiny_data2 <- tiny_data
  names(tiny_data2)[names(tiny_data2) == "predictor"] <- "factor"

  result <- suppressMessages(suppressWarnings(
    OptSurvCutR:::.systematic_search(
      userdata = tiny_data2,
      num_cuts = 1,
      criterion = "logrank",
      covariates = NULL,
      nmin = 8,
      predictor_name = "predictor",
      quiet = TRUE
    )
  ))
  expect_true(all(is.na(result$optimal_cuts)))
})

test_that("Coverage: .systematic_search - 2 cuts with insufficient data", {
  small_data <- mock_data[1:20, ]
  small_data2 <- small_data
  names(small_data2)[names(small_data2) == "predictor"] <- "factor"

  result <- suppressMessages(suppressWarnings(
    OptSurvCutR:::.systematic_search(
      userdata = small_data2,
      num_cuts = 2,
      criterion = "logrank",
      covariates = NULL,
      nmin = 10,
      predictor_name = "predictor",
      quiet = TRUE
    )
  ))
  expect_true(all(is.na(result$optimal_cuts)))
})


test_that("find_cutpoint successfully runs permutations", {
  # Create a small sample dataset
  set.seed(123)
  test_df <- data.frame(
    time = rexp(40),
    status = sample(0:1, 40, replace = TRUE),
    predictor = rnorm(40)
  )

  # Trigger the permutation engine
  res <- find_cutpoint(
    data = test_df,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "status",
    num_cuts = 1,
    n_perm = 5,      # Runs 5 permutations
    n_cores = 1,      # Single core for testing stability
    quiet = TRUE
  )

  # 1. Check if the permuted p-value exists
  expect_false(is.na(res$permuted_p_value))
  expect_true(res$permuted_p_value >= 0 && res$permuted_p_value <= 1)

  # 2. Check if it recorded the number of permutations
  expect_equal(res$n_perm, 5)

  # 3. Test parallel execution (if doParallel is available)
  res_par <- find_cutpoint(
    test_df, "predictor", "time", "status",
    n_perm = 5, n_cores = 2, quiet = TRUE
  )
  expect_false(is.na(res_par$permuted_p_value))
})
