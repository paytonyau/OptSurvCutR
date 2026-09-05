# Find Optimal Number of Cut-points for Survival Data

Finds optimal cut-point number (0 to \`max_cuts\`) for a Cox model by
comparing AIC, AICc, or BIC. Features hardware-accelerated grouping
iterations via Rcpp compilation hooks and robust UX constraint warnings.

## Usage

``` r
find_cutpoint_number(
  data,
  predictor,
  outcome_time,
  outcome_event,
  method = "systematic",
  criterion = "BIC",
  covariates = NULL,
  max_cuts = 2,
  nmin = 0.1,
  seed = NULL,
  max.generations = NULL,
  pop.size = NULL,
  use_cpp = TRUE,
  ...
)

# S3 method for class 'find_cutpoint_number_result'
print(x, ...)

# S3 method for class 'find_cutpoint_number_result'
summary(
  object,
  show_comparison_table = TRUE,
  show_best_model_details = TRUE,
  show_group_counts = TRUE,
  show_medians = TRUE,
  show_ph_test = TRUE,
  plot.it = FALSE,
  ...
)

# S3 method for class 'find_cutpoint_number_result'
plot(x, y, ...)

# S3 method for class 'find_cutpoint_number_result'
print(x, ...)

# S3 method for class 'find_cutpoint_number_result'
summary(
  object,
  show_comparison_table = TRUE,
  show_best_model_details = TRUE,
  show_group_counts = TRUE,
  show_medians = TRUE,
  show_ph_test = TRUE,
  plot.it = FALSE,
  ...
)

# S3 method for class 'find_cutpoint_number_result'
plot(x, y, ...)
```

## Arguments

- data:

  Input data frame.

- predictor:

  Continuous predictor variable name (character).

- outcome_time:

  Time-to-event variable name (character).

- outcome_event:

  Event indicator name (0/1) (character).

- method:

  \`"systematic"\` (max_cuts \<= 2) or \`"genetic"\`.

- criterion:

  \`"AIC"\`, \`"AICc"\` or \`"BIC"\`.

- covariates:

  Character vector of covariate names (optional).

- max_cuts:

  Max number of cut-points to test (non-negative int).

- nmin:

  Min. group size (count or proportion).

- seed:

  Integer or \`NULL\`; random seed for \`rgenoud\`.

- max.generations:

  Integer; generations for \`rgenoud\`. If \`NULL\`, dynamically scales.

- pop.size:

  Integer; population size for \`rgenoud\`. If \`NULL\`, dynamically
  scales.

- use_cpp:

  Logical. Automatically checks and calls compiled C++ routines via
  \`Rcpp\`. Default is \`TRUE\`.

- ...:

  Additional arguments passed down to downstream rendering pipelines.

- x:

  A `find_cutpoint_number_result` object.

- object:

  A `find_cutpoint_number_result` object for analysis overview.

- show_comparison_table:

  Logical. Show information criteria comparison matrix?

- show_best_model_details:

  Logical. Show descriptive layers for the optimal selection?

- show_group_counts:

  Logical. Show categorised patient split breakdowns?

- show_medians:

  Logical. Show Kaplan-Meier time threshold tracking?

- show_ph_test:

  Logical. Display Schoenfeld residuals test?

- plot.it:

  Logical. If `TRUE`, automatically prints the information criterion
  line chart.

- y:

  Unused mandatory base parameter required for graphic dispatcher
  pairing inheritance.

## Value

An S3 object (\`find_cutpoint_number_result\`) with \`results\`,
\`parameters\`, \`userdata\`, \`optimal_num_cuts\`, \`optimal_cuts\`,
and \`candidate_cuts\`.

## Details

\`method = "systematic"\`: grid search respecting \`nmin\`. \`method =
"genetic"\`: \`rgenoud\` global optimisation. Systematic search is slow
for \`max_cuts \> 2\`; use \`genetic\`. Core vector partitions are
calculated in compiled C++ via \`Rcpp\` for optimal performance.

## srrstats compliance

.

.

.

.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)

  # Generate a pristine simulated clinical tracking baseline template
  set.seed(42)
  sim_data <- data.frame(
    time = rexp(40, rate = 0.1),
    event = sample(c(0, 1), 40, replace = TRUE),
    biomarker = rnorm(40, mean = 6, sd = 1.2)
  )

  # Sweep information criteria fit columns up to a 2-cut matrix max
  num_fit <- find_cutpoint_number(
    data = sim_data,
    predictor = "biomarker",
    outcome_time = "time",
    outcome_event = "event",
    max_cuts = 2,
    method = "systematic",
    criterion = "BIC",
    nmin = 5,
    quiet = TRUE
  )
  summary(num_fit)
}
#> ℹ Finding optimal cut number: method = systematic
#> ℹ Profiling IC surface for 1 cut-point(s)...
#> No valid cut-points found for 1 cut(s).
#> ℹ Profiling IC surface for 2 cut-point(s)...
#> No valid cut-points found for 2 cut(s).
#> ! All tested model cut-points violated localised subgroup size constraints during runtime search iterations.
#> 
#> ── Optimal Cut-point Number Analysis (Systematic) ──────────────────────────────
#> ✔ Best Model: 0 Cut-points (Criterion: BIC)
#> 
#> 
#> ── 1. Model Comparison ──
#> 
#>  Marker num_cuts   BIC Delta_BIC BIC_Weight    Evidence
#>       >        0 120.5         0       100% Substantial
#>                1    NA        NA        NA%        <NA>
#>                2    NA        NA        NA%        <NA>
#> 
#> ── 3. Cox Proportional-Hazards ──
#> 
#>   Group    HR Lower Upper P_Value Signif
#>  factor 1.621 1.034 2.543   0.035      *
#> 
#> ℹ Overall Model: Concordance = 0.679 | Log-rank p = 0.035
#> 
#> 
#> ── 4. Time-Dependent Diagnostics (Schoenfeld) ──
#> 
#> ✔ Passed: The proportional hazards assumption holds across the follow-up period (Global p = 0.17).
#> 
#> 
#> ── 5. Analysis Parameters ──
#> 
#> * Search Method: Systematic
#> * Predictor: biomarker
#> * Criterion: BIC
#> * Maximum Cuts: 2
#> * Minimum Group Size (nmin): 5
#> 
```
