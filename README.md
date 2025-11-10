# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml)
  [![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
  [![Codecov](https://codecov.io/gh/paytonyau/OptSurvCutR/branch/main/graph/badge.svg)](https://app.codecov.io/gh/paytonyau/OptSurvCutR/new)
  [![License: GPL-3](https://img.shields.io/badge/License-GPL%203-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.html)
    <!-- badges: end -->
  
`OptSurvCutR` (Optimal Survival Cut-points) is an R package for optimising cut-points in survival analysis, designed for biostatisticians analysing time-to-event data with continuous predictors (e.g., virome abundances in TCGA datasets). It provides a robust workflow to determine the optimal number and location of cut-points, moving beyond median splits to capture non-linear relationships (e.g., U-shaped effects).

## Why OptSurvCutR?
- **Beyond Median Splits**: Identifies the optimal number and location of cut-points using AIC, AICc, or BIC, revealing complex predictor effects.
- **Complete Workflow**: Integrates `find_cutpoint_number()`, `find_cutpoint()`, and `validate_cutpoint()` for end-to-end analysis.
- **Flexible Algorithms**: Offers systematic grid search and genetic algorithms (via `rgenoud`) for efficient multi-cut optimisation.
- **Robust Validation**: Assesses cut-point stability using bootstrap resampling, providing 95% confidence intervals to gauge reliability.
- **User-Friendly**: Provides clear S3 methods (`print`, `summary`, `plot`) for easy interpretation of results.

## Installation
You can install the development version of OptSurvCutR from GitHub. Note that the genetic algorithm (`method = "genetic"`) requires the `rgenoud` package, which should be installed separately from CRAN if you plan to use it.
```r
# Install dependencies
install.packages(c("remotes", "rgenoud", "survival"))

# Install OptSurvCutR
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Quick Workflow with CRC Virome Data
Here is a short example demonstrating the core workflow using the built-in colorectal cancer virome dataset.
```r
# Load necessary packages
library(OptSurvCutR); library(dplyr); library(survival)

# --- 1. Load and prepare the built-in CRC dataset ---
data("crc_virome", package = "OptSurvCutR")

# A quick preparation to make the status column numeric (0=LIVING, 1=DECEASED)
crc_data <- crc_virome %>%
  select(
    time = time_months,
    status_char = status,
    Enterovirus
  ) %>%
  mutate(
    status = as.numeric(substr(status_char, 1, 1))
  ) %>%
  # Remove any rows with missing data
  na.omit()

# --- 2. Find the optimal NUMBER of cut-points ---
# We will test for 0, 1, or 2 cuts using a fast systematic search
number_result <- find_cutpoint_number(
  data = crc_data,
  predictor = "Enterovirus",
  outcome_time = "time",
  outcome_event = "status",
  method = "systematic", # "systematic" is fast for a README
  max_cuts = 2,
  nmin = 0.15, # Ensure groups have at least 15% of subjects
  seed = 42
)

print(number_result)

# The BIC suggests 2 cut-points are optimal for this data.

# --- 3. Find the optimal VALUE of those cut-points ---
# We will find the locations for the 2 optimal cuts
cutpoint_result <- find_cutpoint(
  data = crc_data,
  predictor = "Enterovirus",
  outcome_time = "time",
  outcome_event = "status",
  num_cuts = 2, # Use the result from the step above
  method = "systematic",
  nmin = 0.15,
  seed = 123
)

# --- 4. (Optional) Validate cut-point stability ---
# This step runs a bootstrap and can take a few minutes.
# It is recommended for a full analysis but can be skipped for a quick check.
validation_result <- validate_cutpoint(
  cutpoint_result = cutpoint_result,
  num_replicates = 25, # Use >= 500 for a real analysis
  seed = 456
)

summary(validation_result)

# --- 5. Visualise the Result ---
# The plot reveals three distinct risk groups (Low, Medium, High)
# based on Enterovirus abundance.
# We plot the 'validation_result', which shows the survival curves
# using the original optimal cuts found in step 3.
plot(validation_result, type = "outcome")
```

## Workflow Summary
OptSurvCutR provides a three-step workflow for cut-point analysis:

1.  `find_cutpoint_number()`: Determines the statistically optimal *number* of cut-points using information criteria (AIC, AICc, or BIC).
2.  `find_cutpoint()`: Identifies the precise cut-point *locations* using systematic or genetic algorithms, optimising a chosen survival metric (log-rank, HR, p-value).
3.  `validate_cutpoint()`: Assesses the *stability* of the identified cut-points via bootstrap resampling, providing 95% confidence intervals.
## Resources
- **Vignettes**: See browseVignettes("OptSurvCutR") for detailed tutorials, including analyses of the germination and crc_virome datasets.
- **Package Website**: Full function documentation and articles available at https://paytonyau.github.io/OptSurvCutR/ (or run pkgdown::build_site() locally).
- **Manuscript**: Read the accompanying paper for methodological details and further case studies: Yau, Payton T. O. "OptSurvCutR: Validated Cut-point Selection for Survival Analysis." bioRxiv preprint, posted October 10, 2025. https://doi.org/10.1101/2025.10.08.681246.
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
