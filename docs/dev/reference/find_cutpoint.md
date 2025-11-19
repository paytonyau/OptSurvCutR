# Find Optimal Cut-points for Survival Data

Finds optimal cut-point(s) for a continuous predictor in a time-to-event
(survival) analysis. Uses systematic search (1–2 cuts) or a genetic
algorithm (any number of cuts).

## Usage

``` r
find_cutpoint(
  data,
  predictor,
  outcome_time,
  outcome_event,
  num_cuts = 1,
  method = "systematic",
  criterion = "logrank",
  covariates = NULL,
  nmin = 20,
  seed = NULL,
  maxiter = 100,
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

# S3 method for class 'find_cutpoint'
plot(x, type = "outcome", reference_group = NULL, ...)
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

  The statistic to optimize: \`"logrank"\` (max), \`"hazard_ratio"\`
  (max), or \`"p_value"\` (min). Note: When covariates are provided, the
  \`"logrank"\` criterion is automatically generalized to the Cox score
  test.

- covariates:

  Character vector of covariate names.

- nmin:

  Min. group size (integer count or proportion).

- seed:

  Optional integer seed for \`"genetic"\` method.

- maxiter:

  Number of generations for genetic algorithm (default 100).

- quiet:

  Logical. If \`TRUE\`, suppresses final print.

- ...:

  Additional arguments passed to \`rgenoud\` (e.g., \`popSize\`).

- x:

  An object from \[find_cutpoint()\].

- object:

  An object from \[find_cutpoint()\].

- show_model:

  Logical. Show final Cox model summary?

- show_group_counts:

  Logical. Show N and event counts by group?

- show_medians:

  Logical. Show median survival by group?

- show_ph_test:

  Logical. Show proportional hazards test?

- show_params:

  Logical. Show original function parameters?

- type:

  Plot type: \`"outcome"\`, \`"distribution"\`, or \`"forest"\`.

- reference_group:

  Reference group for forest plot (e.g., \`"G1"\`).

## Value

An object of class \`find_cutpoint\` containing the optimal cut-points,
statistic, and analysis parameters.

## Details

\`method = "systematic"\`: grid search respecting \`nmin\`. \`method =
"genetic"\`: \`rgenoud\` global optimization. Systematic search is slow
for \`num_cuts \> 2\`; use \`genetic\`.

## srrstats compliance

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
Reports\*, 50(3). https://pubmed.ncbi.nlm.nih.gov/5910392/

Mebane Jr, W. R., & Sekhon, J. S. (2011). Genetic Optimization Using
Derivatives: The rgenoud Package for R. \*Journal of Statistical
Software\*, 42, 1–26.
[doi:10.18637/jss.v042.i11](https://doi.org/10.18637/jss.v042.i11)

## Examples

``` r
data(crc_virome)
res <- find_cutpoint(
  data = head(crc_virome, 50),
  predictor = "Alphapapillomavirus",
  outcome_time = "time_months",
  outcome_event = "status",
  num_cuts = 1,
  method = "systematic"
)
#> ℹ Running systematic search...
#> ℹ Testing for 1 cut-point(s)...
#> ✔ Systematic search complete.
#> 
#> ── Optimal Cut-point Analysis for Survival Data (Systematic) ───────────────────
#> • Predictor: Alphapapillomavirus
#> • Criterion: logrank
#> • Optimal Log-Rank Statistic: 2.8814
#> ✔ Recommended Cut-point(s): 3.764
```
