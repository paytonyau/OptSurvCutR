# Find Optimal Number of Cut-points for Survival Data

Finds optimal cut-point number (0 to \`max_cuts\`) for a Cox model by
comparing AIC, AICc, or BIC. Supports systematic search (\`max_cuts \<=
2\`) and genetic algorithm (\`rgenoud\`).

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
  max.generations = 100,
  pop.size = 100,
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

  Integer; generations for \`rgenoud\` (default 100).

- pop.size:

  Integer; population size for \`rgenoud\` (default 100).

- ...:

  Additional arguments passed to \`rgenoud\`.

- x:

  An object from \[find_cutpoint_number()\].

- object:

  An object from \[find_cutpoint_number()\].

- show_comparison_table:

  Logical. Show model comparison table?

- show_best_model_details:

  Logical. Show details for best model?

- show_group_counts:

  Logical. Show group counts for best model?

- show_medians:

  Logical. Show median survival for best model?

- plot.it:

  Logical. Display model selection plot?

- y:

  Unused.

## Value

An S3 object (\`find_cutpoint_number_result\`) with \`results\`,
\`parameters\`, \`userdata\`, \`optimal_num_cuts\`, and
\`optimal_cuts\`.

## Details

\`method = "systematic"\`: grid search respecting \`nmin\`. \`method =
"genetic"\`: \`rgenoud\` global optimisation. Systematic search is slow
for \`max_cuts \> 2\`; use \`genetic\`.

## srrstats compliance

.

## References

Akaike, H. (1974). A new look at the statistical model identification.
\*IEEE Transactions on Automatic Control\*, \*\*19\*\*(6), 716–723.
[doi:10.1109/TAC.1974.1100705](https://doi.org/10.1109/TAC.1974.1100705)

Chang, C., Hsieh, M.-K., Chang, W.-Y., Chiang, A. J., & Chen, J. (2017).
Determining the optimal number and location of cutoff points with
application to data of cervical cancer. \*PLOS ONE\*, 12(4), e0176231.
[doi:10.1371/journal.pone.0176231](https://doi.org/10.1371/journal.pone.0176231)

Chen, Y., Huang, J., He, X., Gao, Y., Mahara, G., Lin, Z., & Zhang, J.
(2019). A novel approach to determine two optimal cut-points of a
continuous predictor with a U-shaped relationship to hazard ratio in
survival data: Simulation and application. \*BMC Medical Research
Methodology\*, 19(1), 96.
[doi:10.1186/s12874-019-0738-4](https://doi.org/10.1186/s12874-019-0738-4)

Schwarz, G. (1978). Estimating the dimension of a model. \*The Annals of
Statistics\*, \*\*6\*\*(2), 461–464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136)

Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model
selection in small samples. \*Biometrika\*, \*\*76\*\*(2), 297–307.
[doi:10.1093/biomet/76.2.297](https://doi.org/10.1093/biomet/76.2.297)

## Examples

``` r
data(crc_virome)
res <- find_cutpoint_number(
  data = head(crc_virome, 50),
  predictor = "Alphapapillomavirus",
  outcome_time = "time_months",
  outcome_event = "status",
  method = "systematic",
  max_cuts = 1
)
#> ℹ nmin 0.1 is a proportion. Min. group size set to 5.
#> ℹ Finding optimal cut number: method = systematic
#> ℹ Testing for 1 cut-point(s)...
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> Warning: Loglik converged before variable  1 ; coefficient may be infinite. 
#> 
#> ── Optimal Cut-point Number Analysis ───────────────────────────────────────────
#> Method: systematic
#> Criterion: BIC
#>  num_cuts   BIC Delta_BIC BIC_Weight    Evidence cuts
#>         0 58.47         4      11.9%    Moderate   NA
#>         1 54.47         0      88.1% Substantial 4.59
#> ✔ Conclusion: 1 cut-point(s) is best based on BIC.
#> Optimal cuts at: 4.59
#> Hint: Use `summary()` for details, `plot()` to visualise.
```
