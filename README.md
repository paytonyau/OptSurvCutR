# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

  [![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml)
  [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![Codecov](https://codecov.io/gh/paytonyau/OptSurvCutR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/paytonyau/OptSurvCutR/new)
  [![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
  
`OptSurvCutR` (Optimal Survival Cut-points) is a validated workflow for finding the optimal number and location of cut-points in survival analysis. It sesigned for statisticians and researchers working with **time-to-event data** and **continuous predictors** (e.g., gene expression, virome abundance, biomarkers). Moves beyond arbitrary median splits to **data-driven, reproducible cut-points**.

## Why OptSurvCutR?
| Feature | Benefit |
|-------|--------|
| **Optimal number of cuts** | Uses AIC, AICc, or BIC to select 0–k cut-points |
| **Flexible search** | Systematic grid or genetic algorithm (`rgenoud`) |
| **Covariate adjustment** | Control for confounders in cut-point discovery |
| **Bootstrap validation** | 95% CI for cut-point stability |
| **Publication-ready plots** | Kaplan-Meier, optimisation curves, diagnostics |

## Installation
You can install the development version of `OptSurvCutR` from GitHub. Note that the genetic algorithm (`method = "genetic"`) requires the `rgenoud` package, which should be installed separately from CRAN if you plan to use it.
```r
# Install dependencies
install.packages(c("remotes", "rgenoud", "survival"))

# Install from GitHub (development version)
if (!require("remotes")) install.packages("remotes")
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Quick Workflow with CRC Virome Data
Here is a short example demonstrating the core workflow using the built-in colorectal cancer virome dataset.
```r
library(OptSurvCutR)
library(dplyr)
library(survival)

# Load built-in TCGA colorectal cancer virome data
data("crc_virome")

# Prepare data: select predictor, time, and event (status is already numeric 0/1)
crc <- crc_virome %>%
  select(time = time_months, status, Enterovirus) %>%
  na.omit()

# Step 1: How many cut-points are optimal?
num_cuts <- find_cutpoint_number(
  data = crc,
  predictor = "Enterovirus",
  outcome_time = "time",
  outcome_event = "status",
  max_cuts = 2,
  nmin = 0.15,
  seed = 42
)
print(num_cuts)
# BIC suggests 2 cut-points

# Step 2: Find the optimal cut-point values
cuts <- find_cutpoint(
  data = crc,
  predictor = "Enterovirus",
  outcome_time = "time",
  outcome_event = "status",
  num_cuts = 2,
  method = "systematic",
  nmin = 0.15,
  seed = 123
)

# Step 3: Validate stability (bootstrap)
val <- validate_cutpoint(cuts, num_replicates = 25, seed = 456)
summary(val)

# Step 4: Visualise the Result
plot(val, type = "outcome")      # Kaplan-Meier
plot(cuts)                       # Optimization curve
plot_diagnostics(cuts)           # Schoenfeld residuals (PH check)
```

## Workflow Summary
OptSurvCutR provides a three-step workflow for cut-point analysis:

1.  `find_cutpoint_number()`: Determines the statistically optimal *number* of cut-points using information criteria (AIC, AICc, or BIC).
2.  `find_cutpoint()`: Identifies the precise cut-point *locations* using systematic or genetic algorithms, optimising a chosen survival metric (log-rank, HR, p-value).
3.  `validate_cutpoint()`: Assesses the *stability* of the identified cut-points via bootstrap resampling, providing 95% confidence intervals.
## Resources
- **Vignettes**: See browseVignettes("OptSurvCutR") for detailed tutorials, including analyses of the germination and crc_virome datasets.
- **Package Website**: Full function documentation and articles available at https://paytonyau.github.io/OptSurvCutR/.
- **Manuscript**: Yau, Payton T. O. "OptSurvCutR: Validated Cut-point Selection for Survival Analysis." bioRxiv preprint, posted October 10, 2025. https://doi.org/10.1101/2025.10.08.681246.
- **NEWS.md**: See NEWS.md file for recent changes and version history.
- **Code of Conduct**: Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this project, you agree to abide by its terms.

## Citation
If you use OptSurvCutR in your research, please cite the accompanying manuscript:
```r
@Article{,
  author = {Payton T. O. Yau},
  title = {OptSurvCutR: Validated Cut-point Selection for Survival Analysis},
  year = {2025},
  doi = {10.1101/2025.10.08.681246},
  publisher = {Cold Spring Harbor Laboratory},
  url = {[https://www.biorxiv.org/content/10.1101/2025.10.08.681246](https://www.biorxiv.org/content/10.1101/2025.10.08.681246)},
  journal = {bioRxiv}
}
```
A JOSS submission is planned post-rOpenSci review.

## Support OptSurvCutR
If you find `OptSurvCutR` helpful in your survival analysis research, please consider supporting its ongoing development/maintenance without any dedicated funding. Your contribution, big or small, directly helps dedicate more time to keeping the project alive and improving.

<a href="https://buymeacoffee.com/payton.yau" target="_blank">
 <img src="https://img.buymeacoffee.com/button-api/?text=Buy%20me%20a%20coffee&emoji=&slug=payton.yau&button_colour=FFDD00&font_colour=000000&font_family=Poppins&outline_colour=000000&coffee_colour=ffffff" alt="Buy Me A Coffee" width="150"></a>


## License
Licensed under the GPL-3 License.

## Contact
For questions or feedback, open an issue at [GitHub Issues](https://github.com/paytonyau/OptSurvCutR/issues).
