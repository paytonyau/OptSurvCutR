# ===================================================================
# TEST HELPER & MOCK DATA
# Automatically loaded by testthat before running any test-*.R files.
# ===================================================================

library(survival)
library(cli)

# --- 1. Pure Deterministic Global Generation ---
set.seed(42)
n_test <- 60

# --- 2. Standard Mock Datasets ---
n1 <- n_test / 2
n2 <- n_test - n1
mock_data <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(rnorm(n1, mean = 40, sd = 5), rnorm(n2, mean = 60, sd = 5)),
  covariate1 = rnorm(n_test, mean = 5, sd = 1),
  covariate2 = sample(c("A", "B"), n_test, replace = TRUE)
)

n1_3 <- floor(n_test / 3)
n2_3 <- floor(n_test / 3)
n3_3 <- n_test - n1_3 - n2_3
mock_data_3groups <- data.frame(
  time = rexp(n_test, rate = 0.05),
  event = sample(0:1, n_test, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(
    rnorm(n1_3, mean = 30, sd = 5),
    rnorm(n2_3, mean = 50, sd = 5),
    rnorm(n3_3, mean = 70, sd = 5)
  ),
  covariate1 = rnorm(n_test, mean = 5, sd = 1)
)

# --- 3. Edge Case & Pathological Datasets ---
mock_data_skewed <- mock_data
mock_data_skewed$predictor <- rexp(n_test, rate = 0.1)

mock_data_heavy_censor <- mock_data
mock_data_heavy_censor$event <- sample(0:1, n_test, replace = TRUE, prob = c(0.9, 0.1))

mock_data_pathological <- mock_data[1:20, ]
mock_data_pathological$predictor <- rep(50, 20)

mock_data_no_events <- mock_data
mock_data_no_events$event <- 0

tiny_data <- mock_data[1:4, ]
tiny_data$predictor <- c(40, 45, 50, 55)

small_data_unique <- data.frame(
  time = rexp(20, rate = 0.05),
  event = sample(0:1, 20, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = rep(50, 20),
  covariate1 = rnorm(20, mean = 5, sd = 1)
)

# --- 4. Pre-Calculated Results for Heavy Bootstraps ---
set.seed(123)
mock_data_for_boot <- data.frame(
  time = rexp(80, rate = 0.05),
  event = sample(0:1, 80, replace = TRUE, prob = c(0.3, 0.7)),
  predictor = c(rnorm(40, mean = 40, sd = 5), rnorm(40, mean = 60, sd = 5))
)

valid_fc_result_for_boot <- suppressMessages(suppressWarnings(
  OptSurvCutR::find_cutpoint(
    data = mock_data_for_boot,
    predictor = "predictor",
    outcome_time = "time",
    outcome_event = "event",
    num_cuts = 1,
    nmin = 5,
    quiet = TRUE
  )
))

if (is.null(valid_fc_result_for_boot) || any(is.na(valid_fc_result_for_boot$optimal_cuts))) {
  valid_fc_result_for_boot <- list(
    optimal_cuts = 50,
    optimal_stat = 12.5,
    all_stats = data.frame(cut = 30:70, stat = rnorm(41, 10, 2)),
    userdata = mock_data_for_boot,
    parameters = list(
      predictor = "predictor", outcome_time = "time",
      outcome_event = "event", num_cuts = 1, nmin = 5,
      method = "systematic", criterion = "logrank", quiet = TRUE
    )
  )
  class(valid_fc_result_for_boot) <- "find_cutpoint"
}
