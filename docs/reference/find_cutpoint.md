# Find Optimal Cut-points for Survival Data

Finds optimal cut-point(s) for a continuous predictor in a time-to-event
(survival) analysis. Uses systematic search (1–2 cuts) or a genetic
algorithm (any number of cuts). Features high-speed integer partitioning
via compiled C++ vector assignments and automated quantile grid
downsampling.

## Usage

``` r
find_cutpoint(
  data,
  predictor,
  outcome_time,
  outcome_event,
  num_cuts = 1,
  method = c("systematic", "genetic"),
  criterion = c("logrank", "hazard_ratio", "p_value"),
  covariates = NULL,
  nmin = 20,
  seed = NULL,
  max.generations = NULL,
  pop.size = NULL,
  n_perm = 0,
  n_cores = 1,
  use_cpp = TRUE,
  grid_by = 0.01,
  quiet = FALSE,
  candidate_cuts = NULL,
  ...
)

# S3 method for class 'find_cutpoint'
print(x, ...)

# S3 method for class 'find_cutpoint'
summary(
  object,
  show_model = TRUE,
  show_group_counts = TRUE,
  show_medians = TRUE,
  show_ph_test = TRUE,
  show_params = TRUE,
  ...
)
```

## Arguments

- data:

  A data frame containing the analysis variables.

- predictor:

  The continuous predictor variable name (character).

- outcome_time:

  The time-to-event variable name (character).

- outcome_event:

  The event status variable name (character, 0 or 1).

- num_cuts:

  The number of cut-points to find. Default is 1.

- method:

  Algorithm search type: \`"systematic"\` or \`"genetic"\`.

- criterion:

  The statistic to optimise: \`"logrank"\`, \`"hazard_ratio"\`, or
  \`"p_value"\`.

- covariates:

  Character vector of covariate names (optional).

- nmin:

  Min. group size (integer count or proportion).

- seed:

  Optional integer seed for reproducible genetic search.

- max.generations:

  Max generations for genetic algorithm. If \`NULL\`, dynamically
  scales.

- pop.size:

  Population size for genetic algorithm. If \`NULL\`, dynamically
  scales.

- n_perm:

  Number of permutations to run for an adjusted p-value. Default is 0.

- n_cores:

  Number of CPU cores for parallel permutations. Default is 1.

- use_cpp:

  Logical. Checks and calls compiled C++ routines via \`Rcpp\`. Default
  is \`TRUE\`.

- grid_by:

  Percentile step increment for systematic grid downsampling. Default is
  0.01.

- quiet:

  Logical. If \`TRUE\`, suppresses operational console alerts.

- candidate_cuts:

  Optional vector of pre-filtered cuts defining a narrow search space.

- ...:

  Additional arguments passed down to downstream rendering pipelines.

- x:

  A find_cutpoint result object.

- object:

  A find_cutpoint result object for summary evaluation.

- show_model:

  Logical. Whether to print the full Cox model summary frame.

- show_group_counts:

  Logical. Whether to show stratified sample split counts.

- show_medians:

  Logical. Whether to display Kaplan-Meier median tracking times.

- show_ph_test:

  Logical. Display the proportional hazards validation check.

- show_params:

  Logical. Print original baseline parameters.

## Value

An object of class \`find_cutpoint\` containing the optimal cut-points,
statistic, and analysis parameters.

## Details

\`method = "systematic"\`: grid search respecting \`nmin\`. Optimised
via internal quantiles. \`method = "genetic"\`: \`rgenoud\` global
optimisation. Systematic search is slow for \`num_cuts \> 2\`; use
\`genetic\`. Core vector partitions are calculated in compiled C++ via
\`Rcpp\` for optimal performance.

## srrstats compliance

.

.

.

## References

Altman, D. G., Lausen, B., Sauerbrei, W., & Schumacher, M. (1994).
Dangers of Using “Optimal” Cutpoints in the Evaluation of Prognostic
Factors. \*JNCI: Journal of the National Cancer Institute\*, 86(11),
829–835.
[doi:10.1093/jnci/86.11.829](https://doi.org/10.1093/jnci/86.11.829)

Cox, D. R. (1972). Regression Models and Life-Tables. \*Journal of the
Royal Statistical Society: Series B (Methodological)\*, 34(2), 187–202.
[doi:10.1111/j.2517-6161.1972.tb00899.x](https://doi.org/10.1111/j.2517-6161.1972.tb00899.x)

Mantel, N. (1966). Evaluation of survival data and two new rank order
statistics arising in its consideration. \*Cancer Chemotherapy
Reports\*, 50(3).

Mebane Jr, W. R., & Sekhon, J. S. (2011). Genetic Optimisation Using
Derivatives: The rgenoud Package for R. \*Journal of Statistical
Software\*, 42, 1–26.
[doi:10.18637/jss.v042.i11](https://doi.org/10.18637/jss.v042.i11)

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)

  # Create a lightweight, reproducible simulation baseline dataset
  set.seed(42)
  sim_data <- data.frame(
    time = rexp(30, rate = 0.1),
    event = sample(c(0, 1), 30, replace = TRUE),
    biomarker = rnorm(30, mean = 5, sd = 1.5)
  )

  # Execute an exhaustive systematic threshold discovery sweep
  fit <- find_cutpoint(
    data = sim_data,
    predictor = "biomarker",
    outcome_time = "time",
    outcome_event = "event",
    num_cuts = 1,
    method = "systematic",
    criterion = "logrank",
    nmin = 5,
    quiet = TRUE
  )
  print(fit)
}
#> 
#> ── Optimal Cut-point Analysis for Survival Data (Systematic) ───────────────────
#> • Predictor: biomarker
#> • Criterion: logrank
#> • Optimal Log-Rank Statistic: 12.4395
#> ✔ Recommended Cut-point(s): 4.827
#> Hint: Use `summary()` for clinical details and Cox regression.
```
