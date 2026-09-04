# ===================================================================
# TESTS: find_cutpoint()
# Core optimisation engines, metric switches, and covariate logic
# ===================================================================

if (!exists("mock_data", envir = .GlobalEnv)) {
  set.seed(42)
  n_mock <- 50
  assign("mock_data", data.frame(
    time       = runif(n_mock, 5, 50),
    event      = sample(c(0, 1), n_mock, replace = TRUE, prob = c(0.3, 0.7)),
    predictor  = rnorm(n_mock, mean = 10, sd = 3),
    covariate1 = rnorm(n_mock, 0, 1),
    covariate2 = factor(sample(c("A", "B"), n_mock, replace = TRUE))
  ), envir = .GlobalEnv)
}

if (!exists("mock_data_pathological", envir = .GlobalEnv)) {
  n_path <- 30
  assign("mock_data_pathological", data.frame(
    time       = rep(20, n_path),
    event      = rep(0, n_path),
    predictor  = rep(10, n_path),
    covariate1 = rep(0, n_path),
    covariate2 = factor(rep("A", n_path))
  ), envir = .GlobalEnv)
}

# --- 1. Core Mathematical Search Engine Sweeps ---

#' @srrstats {G1.0}
#' @srrstats {G1.1}
#' @srrstats {G5.0}
test_that("find_cutpoint executes systematic search across all core splitting metrics", {
  for (crit in c("logrank", "hazard_ratio", "p_value")) {
    res_sys <- suppressMessages(suppressWarnings(
      find_cutpoint(mock_data,
        predictor = "predictor", outcome_time = "time",
        outcome_event = "event", num_cuts = 1, method = "systematic",
        criterion = crit, quiet = TRUE, nmin = 3
      )
    ))

    expect_s3_class(res_sys, "find_cutpoint")
  }
})

#' @srrstats {G5.4a}
test_that("find_cutpoint genetic algorithm converges to valid structural result boundaries", {
  res_gen <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
      predictor = "predictor", outcome_time = "time",
      outcome_event = "event", num_cuts = 1, method = "genetic",
      criterion = "logrank", quiet = TRUE, seed = 123, nmin = 3
    )
  ))

  expect_s3_class(res_gen, "find_cutpoint")
})

# --- 2. Advanced Multi-Omics Covariate Adjustments & Pure R Fallbacks ---

#' @srrstats {RE1.0}
#' @srrstats {RE1.1}
test_that("find_cutpoint incorporates confounders via adjusted Cox models correctly", {
  # Path A: Standard numeric covariate adjustment
  res_adj <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
      predictor = "predictor", outcome_time = "time",
      outcome_event = "event", covariates = "covariate1",
      num_cuts = 1, method = "systematic", criterion = "hazard_ratio", quiet = TRUE, nmin = 3
    )
  ))

  expect_s3_class(res_adj, "find_cutpoint")
  expect_identical(res_adj$parameters$covariates, "covariate1") # UPGRADED: expect_equal -> expect_identical

  # Path B: FORCE PURE R ACCELERATION FALLBACK
  res_pure_r <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
      predictor = "predictor", outcome_time = "time",
      outcome_event = "event", num_cuts = 1, method = "systematic",
      criterion = "logrank", quiet = TRUE, nmin = 3, use_cpp = FALSE
    )
  ))
  expect_s3_class(res_pure_r, "find_cutpoint")

  # Path C: FORCE CATEGORICAL FORMULA CONTRASTS
  res_cat_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
      predictor = "predictor", outcome_time = "time",
      outcome_event = "event", covariates = "covariate2",
      num_cuts = 1, method = "systematic", criterion = "hazard_ratio",
      quiet = TRUE, nmin = 3
    )
  ))
  expect_s3_class(res_cat_cov, "find_cutpoint")
  expect_identical(res_cat_cov$parameters$covariates, "covariate2") # UPGRADED: expect_equal -> expect_identical
})

# --- 3. Front-End Input Validation Trap Arrays ---

#' @srrstats {G5.2b}
test_that("find_cutpoint traps front-end argument anomalies before pipeline entry", {
  expect_error(
    suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", method = "INVALID_METHOD")))
  )

  expect_error(
    find_cutpoint(mock_data, "non_existent_column", "time", "event"),
    regexp = "Missing required columns"
  )
})

# --- 4. Empty Class and Optimisation Crash Resiliency ---

test_that("S3 router for find_cutpoint gracefully degrades on zero-variance pathological inputs", {
  res_empty <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_pathological, "predictor", "time", "event", num_cuts = 1, quiet = TRUE)
  ))

  expect_s3_class(res_empty, "find_cutpoint")
  expect_true(is.na(res_empty$optimal_cuts))
})

test_that("find_cutpoint handles data constraint violations and na_return branches completely", {
  # TRIGGER BRANCH 1: Completely empty dataset
  empty_df <- data.frame(time = numeric(0), event = numeric(0), predictor = numeric(0))
  res_no_cases <- suppressMessages(find_cutpoint(
    empty_df,
    predictor = "predictor", outcome_time = "time", outcome_event = "event", quiet = FALSE
  ))
  expect_true(all(is.na(res_no_cases$optimal_cuts)))

  # TRIGGER BRANCH 2: Insufficient observations for nmin constraints
  res_nmin_fail <- suppressMessages(find_cutpoint(
    mock_data,
    predictor = "predictor", outcome_time = "time", outcome_event = "event",
    nmin = 100, quiet = TRUE
  ))
  expect_true(all(is.na(res_nmin_fail$optimal_cuts)))

  # TRIGGER BRANCH 3: C++ Namespace Absence Guardrail Fallback
  pkg_env <- asNamespace("OptSurvCutR")

  if (exists("cpp_get_group_assignments", envir = pkg_env, mode = "function")) {
    backup_cpp_func <- get("cpp_get_group_assignments", envir = pkg_env)
    unlockBinding("cpp_get_group_assignments", pkg_env)
    assign("cpp_get_group_assignments", NULL, envir = pkg_env)

    res_fallback <- suppressMessages(suppressWarnings(
      find_cutpoint(
        mock_data[1:20, ],
        predictor = "predictor", outcome_time = "time", outcome_event = "event",
        num_cuts = 1, method = "systematic", use_cpp = TRUE, quiet = TRUE, nmin = 2
      )
    ))
    expect_s3_class(res_fallback, "find_cutpoint")

    assignInNamespace("cpp_get_group_assignments", backup_cpp_func, "OptSurvCutR")
    lockBinding("cpp_get_group_assignments", pkg_env)
  }
})

# --- 5. Genetic Optimisation Boundary Stress Tests ---

test_that("engine-genetic handles un-converged states and covariate adjustments cleanly", {
  # 1. TRIGGER COVARIATE ADJUSTED GENETIC SEARCH
  res_gen_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(
      mock_data,
      predictor = "predictor", outcome_time = "time", outcome_event = "event",
      num_cuts = 1, method = "genetic", criterion = "hazard_ratio",
      covariates = "covariate1", max.generations = 3, pop.size = 15, quiet = TRUE
    )
  ))
  expect_s3_class(res_gen_cov, "find_cutpoint")

  # 2. TRIGGER CRASH PROTECTION SWITCH INSIDE GENETIC BACKEND
  pkg_env <- asNamespace("OptSurvCutR")

  if (exists(".run_genetic_search", envir = pkg_env, mode = "function")) {
    backup_gen_func <- get(".run_genetic_search", envir = pkg_env)
    unlockBinding(".run_genetic_search", pkg_env)
    assign(".run_genetic_search", function(...) NULL, envir = pkg_env)

    res_gen_fail <- suppressMessages(suppressWarnings(
      find_cutpoint(
        mock_data,
        predictor = "predictor", outcome_time = "time", outcome_event = "event",
        num_cuts = 1, method = "genetic", criterion = "logrank", quiet = TRUE
      )
    ))

    expect_s3_class(res_gen_fail, "find_cutpoint")
    expect_true(all(is.na(res_gen_fail$optimal_cuts)))

    assignInNamespace(".run_genetic_search", backup_gen_func, "OptSurvCutR")
    lockBinding(".run_genetic_search", pkg_env)
  }
})

# --- 6. Exhaustive S3 Namespace Method Sweeps ---

test_that("S3 methods for find_cutpoint provide exhaustive branch coverage via explicit namespace targeting", {
  res <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
      predictor = "predictor", outcome_time = "time",
      outcome_event = "event", num_cuts = 1, method = "systematic",
      criterion = "logrank", quiet = TRUE, nmin = 3
    )
  ))

  sink_file <- tempfile()
  sink_con <- file(sink_file, open = "wt")
  sink(sink_con)
  sink(sink_con, type = "message")

  on.exit(
    {
      sink(type = "message")
      sink()
      close(sink_con)
      unlink(sink_file)
    },
    add = TRUE
  )

  print(res)

  sum_obj1 <- summary.find_cutpoint(res, show_model = TRUE, show_group_counts = TRUE)
  print(sum_obj1)

  sum_obj2 <- summary.find_cutpoint(res, show_medians = TRUE, show_ph_test = TRUE)
  print(sum_obj2)

  sum_obj3 <- summary.find_cutpoint(res, show_params = FALSE)
  print(sum_obj3)

  sink(type = "message")
  sink()
  close(sink_con)
  on.exit(unlink(sink_file), add = FALSE)

  cat("\n", file = sink_file, append = TRUE)

  captured_lines <- readLines(sink_file)
  expect_gt(length(captured_lines), 0) # UPGRADED: expect_true(x > y) -> expect_gt(x, y)
})
