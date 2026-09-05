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
# \donttest{
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
#> ℹ Running regulared systematic search for 1 cut-point(s)...
#> ✔ Systematic grid optimation complete.
#> ℹ Using random seed 123 for reproducibility.
#> ℹ Bootstrap `nmin` not set. Using 9 (90% of original) to improve stability.
#> ℹ Validating 1 cut(s) from 'systematic' search using 'logrank' over regularised coordinate lattice.
#> ℹ Running 25 replicates sequentially (n_cores = 1).
#> ✔ 25 replicates completed.
#> Cut-point Stability Analysis (Bootstrap)
#> ----------------------------------------
#> Original Optimal Cut-point(s): 3.666 
#> Successful Replicates: 25 / 25 ( 100 %)
#> Failed Replicates: 0 
#> 
#> 95% Confidence Intervals
#> ------------------------
#>       Lower Upper
#> Cut 1 3.118 5.252
#> 
#> Bootstrap Summary Statistics
#> ---------------------------
#>      Cut  Mean    SD Median    Q1    Q3
#> 25% Cut1 3.991 0.657  3.671 3.585 4.522
#> 
#> Hint: Use `summary()` or `plot()` to visualise stability.
#> Cut-point Stability Analysis (Bootstrap)
#> ----------------------------------------
#> Original Optimal Cut-point(s): 3.666 
#> 
#> Bootstrap Distribution Summary
#> -----------------------------
#>      Cut  Mean    SD Median    Q1    Q3
#> 25% Cut1 3.991 0.657  3.671 3.585 4.522
#> 
#> 95% Confidence Intervals
#> ------------------------
#>       Lower Upper
#> Cut 1 3.118 5.252
#> 
#> Validation Parameters
#> ---------------------
#> Replicates Requested: 25 
#> Successful Replicates: 25 / 25 ( 100 %)
#> Failed Replicates: 0 
#> Cores Used: 1 
#> Seed: 123 
#> Minimum Group Size (nmin): 9 
#> Method: systematic 
#> Criterion: logrank 
#> Covariates: None 
#> 
#> 
#> Stability Assessment:
#> ---------------------
#> Maximum CI Width (Relative to 10th-90th Percentile Range): 120.7%
#> ✖ Model Status: UNSTABLE (Tier 4)
#> ! The primary source of instability is Cut 1.
#> ✖ Recommendation: Reduce `num_cuts` or increase `nmin`.
#> 

# }
```
