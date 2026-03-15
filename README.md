# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

<!-- badges: start -->
[![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml)
[![Codecov](https://codecov.io/gh/paytonyau/OptSurvCutR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/paytonyau/OptSurvCutR)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

`OptSurvCutR` (**Opt**imal **Surv**ival **Cut**-points **R**) provides a rigorous, reproducible **three-step workflow** for discovering the optimal number and location of cut-points in time-to-event (survival) data. Designed for continuous predictors (e.g., gene expression, virome abundance, biomarkers), it moves beyond arbitrary median splits to fully **data-driven stratification**.

![OptSurvCutR Graphical Abstract](image/graphical_abstract.jpg)

## Why OptSurvCutR?

| Feature                        | Benefit                                                               |
|-------------------------------|-----------------------------------------------------------------------|
| Optimal number of cuts        | Uses AIC, AICc, or BIC to select 0–k cut-points                       |
| Flexible search               | Systematic grid or genetic algorithm (`rgenoud`)                     |
| Covariate adjustment          | Control for confounders during cut-point discovery                    |
| Bootstrap validation          | 95% confidence intervals for cut-point stability                     |
| Publication-ready plots       | Kaplan–Meier, optimisation curves, Schoenfeld residual diagnostics   |

## Installation
You can install the development version of `OptSurvCutR` from GitHub. Note that the genetic algorithm (`method = "genetic"`) requires the `rgenoud` package, which should be installed separately from CRAN if you plan to use it.
```r
# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for >2 cuts
install.packages("rgenoud")

# Install the development version from GitHub
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Quick Workflow with CRC Virome Data

```r
library(OptSurvCutR)
library(survival)

data("crc_virome")

crc <- crc_virome %>%
  select(time = time_months, status, Enterovirus) %>%
  na.omit()

# Step 1: Determine optimal number of cut-points
num_cuts <- find_cutpoint_number(
  data = crc, predictor = "Enterovirus",
  outcome_time = "time", outcome_event = "status",
  max_cuts = 2, nmin = 0.15, seed = 42
)
print(num_cuts)   # BIC suggests 2 cut-points

# Step 2: Find the precise cut-point locations
cuts <- find_cutpoint(
  data = crc, predictor = "Enterovirus",
  outcome_time = "time", outcome_event = "status",
  num_cuts = 2, method = "systematic", nmin = 0.15
)

# Step 3: Validate stability with bootstrap
val <- validate_cutpoint(cuts, num_replicates = 200, seed = 456)
summary(val)

# Step 4: Visualise results
plot(cuts, type = "outcome")        # Kaplan–Meier curves
plot(cuts, type = "distribution")  # Predictor + cuts
plot_schoenfeld(cuts)               # Proportional hazards check
plot(val)                           # Bootstrap density + 95% CI
```

## Workflow Summary
OptSurvCutR provides a three-step workflow for cut-point analysis:

1.  `find_cutpoint_number()`: – selects the statistically optimal number of cut-points using information criteria
2.  `find_cutpoint()`: locates exact cut-point values (systematic or genetic search)
3.  `validate_cutpoint()`: assesses stability via bootstrapping with 95% confidence intervals
## Resources
- **Vignettes**: `browseVignettes("OptSurvCutR")`
- **Package Website**:  https://paytonyau.github.io/OptSurvCutR/
- **Manuscript**: Yau, Payton T. O. "OptSurvCutR: Validated Cut-point Selection for Survival Analysis." bioRxiv preprint, posted October 10, 2025. https://doi.org/10.1101/2025.10.08.681246.

## Citation

```r
@article{yau2025optsurvcutr,
  author  = {Yau, Payton T. O.},
  title   = {OptSurvCutR: Validated Cut-point Selection for Survival Analysis},
  year    = {2025},
  doi     = {10.1101/2025.10.08.681246},
  publisher = {Cold Spring Harbor Laboratory},
  journal = {bioRxiv},
  url     = {https://www.biorxiv.org/content/10.1101/2025.10.08.681246}
}
```
A JOSS submission is planned post-rOpenSci review.

## Support OptSurvCutR
If `OptSurvCutR` helps your research, consider buying me a coffee — it directly supports ongoing maintenance with no dedicated funding.

<a href="https://buymeacoffee.com/payton.yau" target="_blank">
 <img src="https://img.buymeacoffee.com/button-api/?text=Buy%20me%20a%20coffee&emoji=&slug=payton.yau&button_colour=FFDD00&font_colour=000000&font_family=Poppins&outline_colour=000000&coffee_colour=ffffff" alt="Buy Me A Coffee" width="150"></a>


## License
Licensed under the GPL-3 License.

## Contact
Questions, suggestions, or issues? Please open a ticket:
https://github.com/paytonyau/OptSurvCutR/issues
