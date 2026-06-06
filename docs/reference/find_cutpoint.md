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
  max.generations = 100,
  pop.size = 100,
  n_perm = 0,
  n_cores = 1,
  use_cpp = TRUE,
  grid_by = 0.01,
  quiet = FALSE,
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

  The continuous predictor variable.

- outcome_time:

  The time-to-event variable.

- outcome_event:

  The event status variable (0 or 1).

- num_cuts:

  The number of cut-points to find. Default is 1.

- method:

  Algorithm: \`"systematic"\` or \`"genetic"\`.

- criterion:

  The statistic to optimise: \`"logrank"\` (max), \`"hazard_ratio"\`
  (max), or \`"p_value"\` (min).

- covariates:

  Character vector of covariate names.

- nmin:

  Min. group size (integer count or proportion).

- seed:

  Optional integer seed for reproducible genetic search.

- max.generations:

  Integer; max generations for genetic algorithm (default 100).

- pop.size:

  Integer; population size for genetic algorithm (default 100).

- n_perm:

  Integer. Number of permutations to run for an adjusted p-value.
  Default is 0. Highly recommended for \`num_cuts \>= 2\` to account for
  optimization bias.

- n_cores:

  Integer. Number of CPU cores for parallel permutations. Default is 1.

- use_cpp:

  Logical. Automatically checks and calls compiled C++ routines via
  \`Rcpp\`. Can be overridden if required. Default is \`TRUE\`.

- grid_by:

  Numeric. Percentile step increment for systematic grid downsampling
  (e.g., 0.01 tests every 1st percentile). If \`NULL\`, tests all unique
  values. Default is 0.01.

- quiet:

  Logical. If \`TRUE\`, suppresses final print.

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

Mebane Jr, W. R., & Sekhon, J. S. (2011). Genetic Optimization Using
Derivatives: The rgenoud Package for R. \*Journal of Statistical
Software\*, 42, 1–26.
[doi:10.18637/jss.v042.i11](https://doi.org/10.18637/jss.v042.i11)

## Examples

``` r
if (FALSE) { # \dontrun{
data(crc_virome)
# Fast 1-cut systematic search using downsampled quantile steps
res <- find_cutpoint(
  data = head(crc_virome, 50),
  predictor = "Alphapapillomavirus",
  outcome_time = "time_months",
  outcome_event = "status",
  num_cuts = 1,
  method = "systematic",
  grid_by = 0.01
)
} # }
```
