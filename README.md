# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

[![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: Stable](https://lifecycle.r-lib.org/articles/figures/lifecycle-stable.svg)](https://lifecycle.r-lib.org/articles/stages.html#Stable)
[![Codecov](https://codecov.io/gh/paytonyau/OptSurvCutR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/paytonyau/OptSurvCutR)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
<img src="man/figures/logo.png" align="right" height="165" alt="OptSurvCutR logo" />

`OptSurvCutR` (**Opt**imal **Surv**ival **Cut**-points **R**) provides a rigorous, reproducible, rOpenSci-compliant pipeline for discovering the optimal number and location of stratification cut-points in time-to-event (survival) data. Designed for continuous predictors (e.g., gene expression, microbiome abundance, clinical biomarkers), it moves beyond arbitrary median splits to deliver fully **data-driven, covariate-adjusted stratification**.


## What's New in Version 0.9.9
We have significantly overhauled the validation and diagnostic engines to ensure your discovered thresholds are mathematically stable, rOpenSci-compliant, and publication-ready:
* **Integer Index-Space Mapping:** Replaced continuous floating-point search spaces with a discrete, bounded integer lattice mapped directly to sorted unique data indices. This shifts the engine math from an infinite decimal space to a finite spectrum of actual observations, eliminating micro-decimal overfitting and stochastic seed drift.
* **The 4-Tier Stability Assessment:** Bootstrap validation now automatically grades your cut-points into four tiers (Optimal, Distinct, Caution, Unstable) based on confidence interval width and overlap metrics.
* **Automated Schoenfeld Diagnostics:** The package now flags time-varying effects natively during optimization to ensure your thresholds do not violate the Proportional Hazards assumption.
* **Enhanced Parameter Controls:** Fine-tune the genetic algorithm with `nmin` wedges and soft boundaries to rescue unstable thresholds and prevent overfitting.
* **Continuous 2D Contour Validation Landscapes:** Added native S3 routing support for projecting complex multi-dimensional bootstrap distribution horizons onto smooth contour peaks.


## Why OptSurvCutR?

| Feature                        | Benefit                                                                 |
|--------------------------------|-------------------------------------------------------------------------|
| Optimal number of cuts         | Uses AIC, AICc, or BIC to mathematically select 0–k cut-points          |
| Covariate adjustment           | Prove independent prognostic value by controlling for confounders       |
| 4-Tier Bootstrap validation    | 95% confidence intervals and automated stability grading for thresholds |
| Schoenfeld Diagnostics         | 2-Tier warning system to verify proportional hazards assumptions        |
| Flexible search engine         | Systematic grid or multithreaded Genetic Algorithm (`rgenoud`)          |
| Publication-ready plots        | Kaplan–Meier, distribution curves, forest plots, & 2D topology surfaces |

## Installation
You can install the development version of `OptSurvCutR` directly from GitHub. Note that the genetic algorithm (`method = "genetic"`) requires the `rgenoud` package, which should be installed separately from CRAN if you plan to search for multiple cut-points.

```r
# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for multi-cut genetic optimization
install.packages("rgenoud")

# Install the package from GitHub
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Quick Workflow with Simulated Cohort Data

```r
library(OptSurvCutR)
library(survival)
library(dplyr)

# Generate a reproducible, synthetic survival dataset
set.seed(123)
n <- 200
crc <- tibble(
  abundance = rnorm(n, mean = 5, sd = 2),
  age = rnorm(n, mean = 60, sd = 10),
  # Generate survival times influenced by biomarker abundance
  time = rexp(n, rate = 0.05 * exp(0.3 * (abundance > 5.5) + 0.02 * age)),
  event = sample(c(0, 1), n, replace = TRUE, prob = c(0.3, 0.7))
) %>% filter(time > 0)

# Step 1: Determine the optimal number of cut-points
num_res <- find_cutpoint_number(
  data = crc, predictor = "abundance",
  outcome_time = "time", outcome_event = "event",
  method = "genetic", criterion = "BIC",
  max_cuts = 3, nmin = 0.25, seed = 123
)
summary(num_res)

# Step 2: Find the precise cut-point locations
cut_res <- find_cutpoint(
  data = crc, predictor = "abundance",
  outcome_time = "time", outcome_event = "event",
  method = "genetic", criterion = "logrank",
  num_cuts = num_res$optimal_num_cuts,  # Dynamically pass the result from Step 1
  nmin = 0.275,                         # Fine-tuned for stability
  n_perm = 50, n_cores = 2, seed = 123
)
summary(cut_res)  # Automatically reports Hazard Ratios & Schoenfeld Diagnostics!

# Step 3: Validate stability with bootstrap
val_res <- validate_cutpoint(
  cutpoint_result = cut_res,
  num_replicates = 150, n_cores = 2, seed = 123
)
summary(val_res)  # Automatically grades threshold stability (Tiers 1-4)!

# Step 4: Visualise outcomes via atomic native plots
plot(cut_res, type = "distribution")   # Continuous Predictor Density Split Map
plot(cut_res, type = "outcome")        # Premium Custom Kaplan-Meier Survival Curves
plot(cut_res, type = "forest")         # Hazard Ratio Forest Plot with Cohort Sample Sizes
plot(cut_res, type = "diagnostic")     # Schoenfeld Residual Proportional Hazards Check
plot_validation(val_res, focus_cuts = c(1, 2)) # 2D Contour Elevation Stability Topology
```

## Workflow Summary
`OptSurvCutR` provides a three-step workflow for cut-point analysis:

1. `find_cutpoint_number()`: Selects the statistically optimal number of cut-points using information criteria.
2. `find_cutpoint()`: Locates exact cut-point values (systematic or genetic search) and reports Schoenfeld diagnostics.
3. `validate_cutpoint()`: Assesses threshold stability via bootstrapping and assigns an automated 4-Tier stability grade.

## Resources
- **Vignettes & Tutorials**: Run `browseVignettes("OptSurvCutR")` within your R session to access complete walk-throughs.
- **Troubleshooting & FAQ**: Review the detailed development notes located directly in the `vignettes/troubleshooting.Rmd` source path.
- **Manuscript**: Yau, Payton T. O. "OptSurvCutR: Validated Cut-point Selection for Survival Analysis." bioRxiv preprint, posted October 18, 2025. https://doi.org/10.1101/2025.10.08.681246.

## Citation
```bibtex
@article{yau2025optsurvcutr,
  author    = {Yau, Payton T. O.},
  title     = {OptSurvCutR: Validated Cut-point Selection for Survival Analysis},
  year      = {2025},
  doi       = {10.1101/2025.10.08.681246},
  publisher = {Cold Spring Harbor Laboratory},
  journal   = {bioRxiv},
  url       = {[https://www.biorxiv.org/content/10.1101/2025.10.08.681246](https://www.biorxiv.org/content/10.1101/2025.10.08.681246)}
}
```

## License
Licensed under the GPL-3 License.

## Contact
Questions, suggestions, or issues? Please open a ticket:
https://github.com/paytonyau/OptSurvCutR/issues
```