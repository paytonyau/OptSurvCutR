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

## References

Efron, B. (1979). Bootstrap Methods: Another Look at the Jackknife.
\*The Annals of Statistics\*, 7(1), 1-26.
[doi:10.1214/aos/1176344552](https://doi.org/10.1214/aos/1176344552)

Rota, M., Antolini, L., & Valsecchi, M. G. (2015). Optimal cut-point
definition in biomarkers: The case of censored failure time outcome.
\*BMC Medical Research Methodology\*, 15(1), 24.
[doi:10.1186/s12874-015-0009-y](https://doi.org/10.1186/s12874-015-0009-y)

## Examples

``` r
# Fast validation on small data (runs in < 2 seconds)
data(crc_virome)

fit <- find_cutpoint(
  data = head(crc_virome, 50),
  predictor = "Alphapapillomavirus",
  outcome_time = "time_months",
  outcome_event = "status",
  num_cuts = 1,
  method = "systematic"
)
#> ℹ Running systematic search for 1 cut-point(s)...
#> ✔ Systematic search complete.

if (!any(is.na(fit$optimal_cuts))) {
  val <- validate_cutpoint(fit, num_replicates = 20, seed = 123)
  print(val)
}
#> ℹ Using random seed 123 for reproducibility.
#> ℹ Bootstrap `nmin` not set. Using 18 (90% of original) to improve stability.
#> ℹ Validating 1 cut(s) from 'systematic' search using 'logrank'.
#> ℹ Running 20 replicates sequentially (n_cores = 1).
#> ✔ 20 replicates completed.
#> Cut-point Stability Analysis (Bootstrap)
#> ----------------------------------------
#> Original Optimal Cut-point(s): 3.764 
#> Successful Replicates: 20 / 20 ( 100 %)
#> Failed Replicates: 0 
#> 
#> 95% Confidence Intervals
#> ------------------------
#>       Lower Upper
#> Cut 1 0.826 4.594
#> 
#> Bootstrap Summary Statistics
#> ---------------------------
#>      Cut  Mean    SD Median    Q1    Q3
#> 25% Cut1 2.794 1.485   3.52 1.353 3.764
#> 
#> Hint: Use `summary()` or `plot()` to visualise stability.
```
