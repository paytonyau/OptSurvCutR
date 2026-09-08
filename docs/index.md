# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

[![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle:
Stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#Stable)
[![Codecov](https://codecov.io/gh/paytonyau/OptSurvCutR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/paytonyau/OptSurvCutR)
[![License:
GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
![OptSurvCutR logo](reference/figures/logo.png)

`OptSurvCutR` (**Opt**imal **Surv**ival **Cut**-points **R**) provides a
reproducible framework for identifying the number and location of
patient stratification thresholds in time-to-event (survival) data, and
for assessing whether those thresholds are reproducible. Designed for
continuous predictors (such as gene expression measurements, microbiome
abundance, or clinical biomarkers), this package moves beyond arbitrary
median splits to deliver **data-driven, covariate-adjusted
stratification**. It is developed against rOpenSci software standards
and is currently under rOpenSci review.

## What’s New in Version 0.11

**Important bug fix.** In versions 0.10.1 and earlier,
[`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md)
did not carry the predictor into its result, so
[`summary()`](https://rdrr.io/r/base/summary.html) computed relative
confidence interval widths against the bootstrap distribution of the
first cut-point rather than the predictor’s 10th–90th percentile range.
Reported widths were inflated and **Stability Tier 1 was unreachable for
any input**. If you have reported a stability tier from an earlier
version, please recompute it.

- **Corrected stability metric:** Relative interval width is now
  computed against the predictor’s percentile spread. Where the
  predictor is unavailable, the width is reported as undefined rather
  than substituted with an unrelated quantity.

- **Machine-readable stability output:**
  [`summary()`](https://rdrr.io/r/base/summary.html) attaches a
  `stability` list to the returned object (`tier`, `percent`,
  `relative_widths` per cut-point, `data_spread`, `separated`,
  `worst_cut`), so tiers can be read programmatically instead of parsed
  from console text.

- **Single cut-point grading:** Interval separation is undefined with
  one boundary, so single-threshold models are now graded on relative
  width alone and can no longer be assigned Tier 3.

- **Integer Index-Space Mapping:** The search operates over a discrete,
  bounded integer lattice mapped to sorted unique observations rather
  than a continuous floating-point space, restricting candidate
  thresholds to observed values.

- **The Four-Tier Stability Assessment:** Bootstrap validation grades
  thresholds into four tiers (Optimal, Distinct, Caution, Unstable) from
  confidence interval width and overlap. The tiers summarise two
  independent diagnostics rather than forming an ordinal scale, and the
  30% and 60% boundaries are pragmatic conventions rather than
  calibrated thresholds.

- **Automated Schoenfeld Diagnostics:**
  [`summary()`](https://rdrr.io/r/base/summary.html) reports a test of
  the proportional hazards assumption alongside the fitted model.

- **2D Contour Validation Landscapes:**
  [`plot_validation()`](https://paytonyau.github.io/OptSurvCutR/reference/plot_validation.md)
  projects the joint bootstrap distribution of two thresholds onto a
  contour map, showing whether they vary independently or together.

## Why OptSurvCutR?

| Feature | Benefit |
|:---|:---|
| **Optimal number of cuts** | Uses AIC, AICc, or BIC to select between $`0`$ and $`k`$ cut-points. |
| **Covariate adjustment** | Selects thresholds conditionally on clinical confounders rather than marginally. |
| **Four-tier bootstrap validation** | Generates 95% confidence intervals for each threshold and an automated stability grade. |
| **Schoenfeld diagnostics** | Reports a two-tier test of the proportional hazards assumption for the fitted model. |
| **Flexible search engines** | Supports a deterministic systematic grid or a multithreaded genetic algorithm (`rgenoud`). |
| **Publication-ready plots** | Renders Kaplan–Meier curves, distribution splits, forest plots, objective surfaces, and 2D stability maps. |

## Installation

Install the released version of `OptSurvCutR` from CRAN:

``` r

install.packages("OptSurvCutR")
```

Alternatively, install the development version from GitHub:

``` r

if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Serum Bilirubin in Primary Biliary Cholangitis

This example uses the Mayo Clinic PBC trial cohort supplied with the
`survival` package, and reproduces the analysis reported in the
manuscript.

``` r

library(OptSurvCutR)
library(survival)

data("pbc", package = "survival")
pbc_clean <- na.omit(pbc[, c("time", "status", "bili", "age", "sex", "edema")])

# Composite endpoint: transplant-free survival
# status: 0 = censored, 1 = transplant, 2 = death
pbc_clean$event <- as.integer(pbc_clean$status %in% c(1, 2))

# Step 1: How many cut-points does the data support?
num_res <- find_cutpoint_number(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "genetic", criterion = "AIC",
  max_cuts = 5, nmin = 0.15, seed = 123
)
summary(num_res)   # AIC minimised at 3 cut-points

# Step 2: Where are the boundaries, adjusting for covariates?
cut_res <- find_cutpoint(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "genetic", criterion = "logrank",
  num_cuts = num_res$optimal_num_cuts,   # carried forward from Step 1
  nmin = 0.15,
  n_perm = 200,        # use >= 1000 when reporting
  n_cores = 2, seed = 123
)
summary(cut_res)   # thresholds, adjusted hazard ratios, Schoenfeld test

# Step 3: Are those boundaries reproducible?
val_res <- validate_cutpoint(
  cutpoint_result = cut_res,
  num_replicates = 200,   # use >= 500 when reporting
  n_cores = 2, seed = 123
)
summary(val_res)   # Tier 3: adjacent confidence intervals overlap

# The stability result is also available programmatically
summary(val_res)$stability$tier

# Step 4: Act on it. A single threshold proves highly reproducible.
cut_1 <- find_cutpoint(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "systematic", criterion = "logrank",
  num_cuts = 1, nmin = 0.15, n_perm = 200,
  n_cores = 2, seed = 123
)
summary(validate_cutpoint(cut_1, num_replicates = 200,
                          n_cores = 2, seed = 123))
# Tier 1: threshold 2.3 mg/dL, bootstrap interval 2.0-3.4

# Visualisation via native S3 plot routes
plot(cut_res, type = "distribution")            # predictor density with thresholds
plot(cut_res, type = "outcome")                 # adjusted Kaplan-Meier curves
plot(cut_res, type = "forest")                  # adjusted hazard ratios
plot(cut_res, type = "diagnostic")              # Schoenfeld residual check
plot_validation(val_res, focus_cuts = c(1, 2))  # joint bootstrap distribution
```

**What this example shows.** Information criteria selected four risk
groups with 98.4% of the AIC weight and adjusted hazard ratios up to
20.9, a result that would pass any test of model fit or statistical
significance. Bootstrap validation showed the boundaries were not
reproducible across resampled cohorts. A single threshold at 2.3 mg/dL
was. Statistical significance and threshold stability are distinct
properties, and Step 3 is what separates them.

See
[`vignette("bilirubin", package = "OptSurvCutR")`](https://paytonyau.github.io/OptSurvCutR/articles/bilirubin.md)
for the full analysis.

## Workflow Summary

`OptSurvCutR` establishes a structured, three-step workflow for
cut-point analysis:

1.  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md):
    Identifies the statistically optimal number of thresholds using
    information criteria.
2.  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md):
    Localises exact cut-point coordinates via systematic or genetic
    search, and reports Schoenfeld diagnostics.
3.  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md):
    Evaluates threshold stability using bootstrap resampling and assigns
    an automated four-tier stability grade.

## Contributing

Contributions, bug reports, and pull requests are welcome! Please read
our [Contributing
Guide](https://paytonyau.github.io/OptSurvCutR/CONTRIBUTING.md) for
details on code standards, local test workflows, and issue reporting.

## Resources

- **Vignettes & Tutorials**: Run `browseVignettes("OptSurvCutR")` within
  your R session to access complete walk-throughs.
- **Manuscript**: Yau, Payton T. O. “OptSurvCutR: Validated Cut-point
  Selection for Survival Analysis.” bioRxiv preprint, posted October
  18, 2025. <https://doi.org/10.1101/2025.10.08.681246>.

## Citation

``` bibtex
@article{yau2025optsurvcutr,
  author    = {Yau, Payton T. O.},
  title     = {OptSurvCutR: Validated Cut-point Selection for Survival Analysis},
  year      = {2025},
  doi       = {10.1101/2025.10.08.681246},
  publisher = {Cold Spring Harbor Laboratory},
  journal   = {bioRxiv},
  url       = {https://www.biorxiv.org/content/10.1101/2025.10.08.681246}
}
```

## License

Licensed under the GPL-3 License.

## Contact

For questions, feature suggestions, or bug reports, please open an issue
tracking ticket: <https://github.com/paytonyau/OptSurvCutR/issues>
