# Validate an Optimal Cut-point Using Bootstrapping

Assesses cut-point stability from
[`find_cutpoint`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
via bootstrap analysis, generating 95% confidence intervals. Streamlined
for survival (time-to-event) analysis.

## Usage

``` r
validate_cutpoint(
  cutpoint_result,
  num_replicates = 500,
  n_cores = 1,
  seed = NULL,
  nmin = NULL,
  ...
)
```

## Arguments

- cutpoint_result:

  An object from
  [`find_cutpoint`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md).

- num_replicates:

  Number of bootstrap replicates. Default is 500.

- n_cores:

  Number of CPU cores to use. Default is 1 (sequential). Set to \> 1 to
  enable parallel processing.

- seed:

  Optional integer for reproducible results.

- nmin:

  Minimum group size for bootstrap runs. Defaults to 90% of original
  `nmin` to reduce failures.

- ...:

  Additional arguments passed to
  [`find_cutpoint`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md)
  (e.g., `pop.size`, `max.generations` for genetic algorithm).

## Value

An object of class `validate_cutpoint_result` with original cuts, 95%
CIs, bootstrap distribution, and parameters.

## srrstats compliance

.

## Examples

``` r
if (FALSE) { # \dontrun{
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)

  # 1. Create a tiny simulated baseline clinical cohort dataset
  set.seed(123)
  n <- 45
  toy_data <- data.frame(
    time = rexp(n, rate = 0.05),
    event = sample(c(0, 1), n, replace = TRUE, prob = c(0.4, 0.6)),
    marker = rnorm(n, mean = 4, sd = 1.2)
  )

  # 2. Locate initial baseline cut-points via systematic search
  initial_cut <- find_cutpoint(
    data = toy_data, predictor = "marker",
    outcome_time = "time", outcome_event = "event",
    num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 10
  )

  # 3. Run a lightweight bootstrap validation stress-test execution loop
  val_res <- validate_cutpoint(
    cutpoint_result = initial_cut,
    num_replicates = 25, # Small iteration tier for rapid check verification
    n_cores = 1,
    seed = 123
  )

  # 4. Invoke structural S3 output verification hooks
  print(val_res)
  summary(val_res)
  plot(val_res)
}
} # }
```
