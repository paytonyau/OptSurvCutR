# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

[![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: Stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#Stable)
[![Codecov](https://codecov.io/gh/paytonyau/OptSurvCutR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/paytonyau/OptSurvCutR)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
<img src="man/figures/logo.png" align="right" height="165" alt="OptSurvCutR logo" />

`OptSurvCutR` finds data-driven, covariate-adjusted thresholds for stratifying patients in survival (time-to-event) data — for biomarkers, gene expression, or any continuous predictor — and tells you whether those thresholds would hold up in a different sample. It is developed against rOpenSci software standards and is currently under [rOpenSci review](https://github.com/ropensci/software-review).

**Most threshold-selection tools stop once a split is statistically significant. OptSurvCutR's third step is what happens next: does the same threshold reappear under resampling, or was it a feature of this one dataset?**

## Why OptSurvCutR?

| Feature | Benefit |
| :--- | :--- |
| **Optimal number of cuts** | AIC, AICc, or BIC choose between 0 and *k* cut-points. |
| **Covariate adjustment** | Thresholds are selected conditional on clinical confounders, not marginally. |
| **Bootstrap stability grade** | Every threshold gets a 95% CI and an automated Tier 1–4 grade. |
| **Schoenfeld diagnostics** | Proportional hazards is checked automatically alongside the fit. |
| **Two search engines** | A deterministic grid, or a multithreaded genetic search (`rgenoud`). |
| **Publication-ready plots** | KM curves, distribution splits, forest plots, objective surfaces, 2D stability maps. |

## Installation

```r
install.packages("OptSurvCutR")

# development version
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Serum Bilirubin in Primary Biliary Cholangitis

Using the Mayo Clinic PBC trial cohort shipped with `survival`, reproducing the manuscript's analysis:

```r
library(OptSurvCutR)
library(survival)

data("pbc", package = "survival")
pbc_clean <- na.omit(pbc[, c("time", "status", "bili", "age", "sex", "edema")])

# Composite endpoint: transplant-free survival (0 = censored, 1 = transplant, 2 = death)
pbc_clean$event <- as.integer(pbc_clean$status %in% c(1, 2))

# 1. How many cut-points does the data support?
num_res <- find_cutpoint_number(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "genetic", criterion = "AIC",
  max_cuts = 5, nmin = 0.15, seed = 123
)
summary(num_res)   # AIC minimised at 3 cut-points

# 2. Where are the boundaries, adjusting for covariates?
cut_res <- find_cutpoint(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "genetic", criterion = "logrank",
  num_cuts = num_res$optimal_num_cuts,   # carried forward, not re-chosen
  nmin = 0.15, n_perm = 1000, n_cores = 2, seed = 123
)
summary(cut_res)   # thresholds, adjusted hazard ratios, Schoenfeld test

# 3. Are those boundaries reproducible?
val_res <- validate_cutpoint(
  cutpoint_result = cut_res,
  num_replicates = 500, n_cores = 2, seed = 123
)
summary(val_res)                       # Tier 3 (OVERLAPPING): adjacent intervals overlap
summary(val_res)$stability$tier        # same result, programmatically

# 4. Act on it - a single threshold proves highly reproducible
cut_1 <- find_cutpoint(
  data = pbc_clean, predictor = "bili",
  outcome_time = "time", outcome_event = "event",
  covariates = c("age", "sex", "edema"),
  method = "systematic", criterion = "logrank",
  num_cuts = 1, nmin = 0.15, n_perm = 1000, n_cores = 2, seed = 123
)
summary(validate_cutpoint(cut_1, num_replicates = 500, n_cores = 2, seed = 123))
# Tier 1: threshold 2.3 mg/dL, bootstrap interval 2.0-3.4

# Visualisation
plot(cut_res, type = "distribution")            # predictor density with thresholds
plot(cut_res, type = "outcome")                 # adjusted Kaplan-Meier curves
plot(cut_res, type = "forest")                  # adjusted hazard ratios
plot(cut_res, type = "diagnostic")               # Schoenfeld residual check
plot_validation(val_res, focus_cuts = c(1, 2))   # joint bootstrap distribution
```

**What this shows.** The 4-group model has 98.4% of the AIC weight and hazard ratios up to 20.9 - it would pass any significance test. Bootstrap validation showed its boundaries don't survive resampling. A single threshold at 2.3 mg/dL does. Significance and stability are different properties; step 3 is what tells them apart.

Full analysis: `vignette("bilirubin", package = "OptSurvCutR")`.

## Workflow

1. **`find_cutpoint_number()`** - how many thresholds does the data support (AIC/AICc/BIC)?
2. **`find_cutpoint()`** - where are they, adjusted for covariates? Includes Schoenfeld diagnostics.
3. **`validate_cutpoint()`** - do they reproduce under bootstrap resampling? Returns a Tier 1-4 grade.

## What's New

### v0.11.1 - two correctness fixes, one label change, one new argument

If you used either combination below under v0.11.0, re-run with v0.11.1.

* **No-covariate model selection:** `find_cutpoint_number(method = "systematic")` could silently return a 0-cut model regardless of the data when `covariates = NULL`, due to a null-model indexing error. Models with covariates were unaffected.
* **Permutation p-value direction:** `criterion = "p_value"` with `n_perm > 0` computed the permutation-adjusted p-value in the wrong direction, understating significance. `"logrank"` and `"hazard_ratio"` were unaffected.
* Also fixed: multi-cut Schoenfeld diagnostic panels collapsing into one mislabelled "Cohort G" facet, and a Kaplan-Meier y-axis label that assumed an overall-survival endpoint.
* **Breaking:** Tier labels `"CAUTION"` and `"DISTINCT"` are renamed to `"OVERLAPPING"` and `"CONSISTENT"`. The numeric `tier` field and tier thresholds are unchanged - this affects only the `tier_label` string and console wording. Code matching on the old label strings needs updating.
* **New:** `validate_cutpoint(..., quiet = TRUE)` suppresses the "Bootstrapping" progress bar - useful in scripts, test suites, or anywhere a live progress bar isn't wanted. Default behaviour is unchanged.

Full technical detail in `NEWS.md`.

### v0.11.0 - stability metric fix

**Recompute any Tier reported from v0.10.1 or earlier.** `validate_cutpoint()` didn't carry the predictor into its result, so relative CI width was computed against the wrong quantity - widths were inflated and **Tier 1 was unreachable for any input**.

* Relative interval width now computed against the predictor's 10th-90th percentile range.
* `summary()` attaches a machine-readable `stability` list (`tier`, `percent`, `relative_widths`, `separated`, ...).
* Single-cut models are graded on width alone (overlap is undefined with one boundary) and can no longer land in Tier 3.
* Search now operates over a discrete lattice of observed values rather than continuous floating point.
* Added automated Schoenfeld diagnostics and 2D contour plots (`plot_validation()`) for joint threshold stability.

## Contributing

Bug reports and pull requests are welcome - see [CONTRIBUTING.md](https://github.com/paytonyau/OptSurvCutR/blob/main/CONTRIBUTING.md) for code standards and local test workflows.

## Resources

- Vignettes: `browseVignettes("OptSurvCutR")`
- Manuscript: Yau, Payton T. O. "OptSurvCutR: Validated Cut-point Selection for Survival Analysis." *bioRxiv*, posted October 18, 2025. https://doi.org/10.1101/2025.10.08.681246

## Citation

```bibtex
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

GPL-3.

## Contact

Questions, feature requests, or bug reports: [open an issue](https://github.com/paytonyau/OptSurvCutR/issues).
