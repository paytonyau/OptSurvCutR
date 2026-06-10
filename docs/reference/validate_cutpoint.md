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
