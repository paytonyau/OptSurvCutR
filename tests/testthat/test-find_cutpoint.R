# ===================================================================
# TESTS: find_cutpoint()
# Core optimization engines, metric switches, and covariate logic
# ===================================================================

# --- 0. LOCAL INFRASTRUCTURE FIXTURE ASSURANCE ---
# Ensures tests execute successfully in isolated covr background threads

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
    event      = rep(0, n_path), # Zero events forces mathematical singularities
    predictor  = rep(10, n_path), # Zero variance predictor
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
      ) # Tight nmin prevents boundary drops
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
  expect_equal(res_adj$parameters$covariates, "covariate1")

  # Path B: FORCE PURE R ACCELERATION FALLBACK
  # Disabling C++ compilation forces coverage on your native R matrix loops
  res_pure_r <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
                  predictor = "predictor", outcome_time = "time",
                  outcome_event = "event", num_cuts = 1, method = "systematic",
                  criterion = "logrank", quiet = TRUE, nmin = 3, use_cpp = FALSE
    )
  ))
  expect_s3_class(res_pure_r, "find_cutpoint")

  # Path C: FORCE CATEGORICAL FORMULA CONTRASTS
  # Passing 'covariate2' (character/factor "A"/"B") exercises your internal formula-building engine
  res_cat_cov <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data,
                  predictor = "predictor", outcome_time = "time",
                  outcome_event = "event", covariates = "covariate2",
                  num_cuts = 1, method = "systematic", criterion = "hazard_ratio",
                  quiet = TRUE, nmin = 3
    )
  ))
  expect_s3_class(res_cat_cov, "find_cutpoint")
  expect_equal(res_cat_cov$parameters$covariates, "covariate2")
})

# --- 3. Front-End Input Validation Trap Arrays ---

#' @srrstats {G5.2b}
test_that("find_cutpoint traps front-end argument anomalies before pipeline entry", {
  # Cleanly intercept base match.arg tracking faults safely
  expect_error(
    suppressMessages(suppressWarnings(find_cutpoint(mock_data, "predictor", "time", "event", method = "INVALID_METHOD")))
  )

  expect_error(
    find_cutpoint(mock_data, "non_existent_column", "time", "event"),
    regexp = "Missing required columns"
  )
})

# --- 4. Empty Class and Optimization Crash Resiliency ---

test_that("S3 router for find_cutpoint gracefully degrades on zero-variance pathological inputs", {
  # Leverage your pre-made global pathological dataset directly
  res_empty <- suppressMessages(suppressWarnings(
    find_cutpoint(mock_data_pathological, "predictor", "time", "event", num_cuts = 1, quiet = TRUE)
  ))

  expect_s3_class(res_empty, "find_cutpoint")
  expect_true(is.na(res_empty$optimal_cuts))
})

test_that("find_cutpoint handles data constraint violations and na_return branches completely", {
  # TRIGGER BRANCH 1: Completely empty dataset (Forces nrow(userdata) == 0 path)
  empty_df <- data.frame(time = numeric(0), event = numeric(0), predictor = numeric(0))
  res_no_cases <- suppressMessages(find_cutpoint(
    empty_df,
    predictor = "predictor", outcome_time = "time", outcome_event = "event", quiet = FALSE
  ))
  expect_true(all(is.na(res_no_cases$optimal_cuts)))

  # TRIGGER BRANCH 2: Insufficient observations for nmin constraints
  # Passing nmin = 100 on a dataset of 50 rows triggers .validate_data_conditions failure
  res_nmin_fail <- suppressMessages(find_cutpoint(
    mock_data,
    predictor = "predictor", outcome_time = "time", outcome_event = "event",
    nmin = 100, quiet = TRUE
  ))
  expect_true(all(is.na(res_nmin_fail$optimal_cuts)))

  # TRIGGER BRANCH 3: C++ Namespace Absence Guardrail Fallback Execution Assurance
  pkg_env <- asNamespace("OptSurvCutR")

  if (exists("cpp_get_group_assignments", envir = pkg_env, mode = "function")) {
    backup_cpp_func <- get("cpp_get_group_assignments", envir = pkg_env)

    # Unlock the package environment namespace to allow mocking during test evaluation
    unlockBinding("cpp_get_group_assignments", pkg_env)

    # Assign NULL so 'exists(mode = "function")' returns FALSE natively
    assign("cpp_get_group_assignments", NULL, envir = pkg_env)

    # Assert that function falls back gracefully and constructs a valid S3 return object
    res_fallback <- suppressMessages(suppressWarnings(
      find_cutpoint(
        mock_data[1:20, ],
        predictor = "predictor", outcome_time = "time", outcome_event = "event",
        num_cuts = 1, method = "systematic", use_cpp = TRUE, quiet = TRUE, nmin = 2
      )
    ))
    expect_s3_class(res_fallback, "find_cutpoint")

    # Restore the Rcpp mapping profile immediately to preserve project stability
    assignInNamespace("cpp_get_group_assignments", backup_cpp_func, "OptSurvCutR")
    lockBinding("cpp_get_group_assignments", pkg_env)
  }
})

# --- 5. Genetic Optimization Boundary Stress Tests ---

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
  # Mock the unexported genetic search engine within the package environment namespace.
  # Forcing it to return NULL mimics a catastrophic optimization or convergence failure!
  pkg_env <- asNamespace("OptSurvCutR")

  if (exists(".run_genetic_search", envir = pkg_env, mode = "function")) {
    backup_gen_func <- get(".run_genetic_search", envir = pkg_env)

    # Unlock the package namespace to inject the mock failure
    unlockBinding(".run_genetic_search", pkg_env)

    # Assign a mock function that returns NULL to simulate a hard crash
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

    # Restore the package environment instantly to maintain suite stability
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

  # Create a temporary null sink connection to catch all console output streams
  sink_file <- tempfile()
  sink_con <- file(sink_file, open = "wt")
  sink(sink_con)
  sink(sink_con, type = "message")

  # Ensure the sink closes cleanly even if assertions throw errors
  on.exit(
    {
      sink(type = "message")
      sink()
      close(sink_con)
      unlink(sink_file)
    },
    add = TRUE
  )

  # 1. Test standard object print routing
  print(res)

  # 2. Test summary object parsing toggles naturally
  sum_obj1 <- OptSurvCutR:::summary.find_cutpoint(res, show_model = TRUE, show_group_counts = TRUE)
  print(sum_obj1)

  sum_obj2 <- OptSurvCutR:::summary.find_cutpoint(res, show_medians = TRUE, show_ph_test = TRUE)
  print(sum_obj2)

  sum_obj3 <- OptSurvCutR:::summary.find_cutpoint(res, show_params = FALSE)
  print(sum_obj3)

  # Close the sinks completely before evaluating results
  sink(type = "message")
  sink()
  close(sink_con)
  on.exit(unlink(sink_file), add = FALSE)

  # Push a trailing token line to clear any EOL warning hooks
  cat("\n", file = sink_file, append = TRUE)

  captured_lines <- readLines(sink_file)
  expect_true(length(captured_lines) > 0)
})
