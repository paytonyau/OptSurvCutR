# OptCutR: An R Package for Optimal Cut-point Discovery and Validation

<!-- Badges will go here once you set up GitHub Actions and submit to CRAN -->
<!-- e.g., [![CRAN status](https://www.r-pkg.org/badges/version/OptCutR)](https://cran.r-project.org/package=OptCutR) -->
<!-- e.g., [![R-CMD-check](https://github.com/your-username/OptCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/your-username/OptCutR/actions/workflows/R-CMD-check.yaml) -->

`OptCutR` provides a comprehensive and flexible workflow to determine, find, and validate optimal cut-points for continuous predictors. It is designed to work with both time-to-event (survival) and binary classification outcomes, making it a versatile tool for researchers in bioinformatics, ecology, medicine, and beyond.

![OptCutR Workflow Diagram](https://i.imgur.com/8J3b2gH.png)

### Key Features

* 🔎 **Complete Discovery Workflow:** A suite of functions guides you from determining the *best number* of cuts, to finding their *exact location*, and finally to *validating* the stability of the result.
* 💪 **Two Powerful Search Methods:**
    * **Systematic Search:** A robust, "brute-force" method that guarantees finding the single best cut-point for both survival and binary outcomes.
    * **Genetic Algorithm:** An efficient optimization algorithm for finding multiple cut-points simultaneously in survival analysis, allowing for stratification into three or more risk groups.
* ✨ **User-Friendly Outputs:** All core functions return objects with custom `print()`, `summary()`, and `plot()` methods, providing clean summaries and publication-ready visualizations with single commands.
* ✅ **Built-in Validation:** Includes a dedicated function for bootstrap validation to assess the stability of discovered cut-points and generate 95% confidence intervals.

## Package Workflow

The `OptCutR` workflow is built around three main functions designed to be used in sequence for a complete analysis:

1.  **`find_cutpoint_number()`**: The first step in an exploratory survival analysis. Use this to get statistical evidence (using AIC or BIC) for the most plausible number of cut-points for your data.
2.  **`find_cutpoint()`**: The main workhorse function. After you know how many cuts to look for (or if you have a specific hypothesis), use this to find their exact locations.
3.  **`validate_cutpoint()`**: The final step. Use this to assess the stability of your discovered cut-point(s) by running a bootstrap analysis and generating confidence intervals.

The table below outlines the recommended functions and settings for each scenario.

| Step | Binary Classification (1 Cut-point) | Survival Analysis (1 Cut-point) | Survival Analysis (2+ Cut-points) |
| :--- | :--- | :--- | :--- |
| **1. Find Number of Cuts** | **Not Applicable.**<br>The decision for 1 cut is hypothesis-driven. | **Optional & Exploratory.**<br>Use `find_cutpoint_number()` to see if the data supports 1 cut. | **Recommended.**<br>Use `find_cutpoint_number()` to determine the most plausible number of cuts. |
| **2. Find Cut-point(s)** | **Recommended Method.**<br>`find_cutpoint(method="systematic")` with majority vote from ROC metrics. | **Recommended Method.**<br>`find_cutpoint(method="systematic")` for the most robust single cut-point. | **Required Method.**<br>`find_cutpoint(method="genetic")` with `num_cuts` set from Step 1. |
| **3. Validate Stability** | Assess stability with bootstrapping.<br>`validate_cutpoint()` | Assess stability with bootstrapping.<br>`validate_cutpoint()` | Assess stability with bootstrapping.<br>`validate_cutpoint()` |

## Installation

You can install the development version of `OptCutR` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("paytonyau/OptCutR")
```

## Quick Example: A Complete Workflow in a Few Lines

Here is a minimal example showing the full workflow for a survival analysis with two cut-points.

```r
# Load the package and some example data
library(OptCutR)
library(survival)

# Step 1: Find the optimal number of cuts (let's assume it's 2)
# num_res <- find_cutpoint_number(data = lung, predictor = "age", outcome_time = "time", outcome_event = "status", max_cuts = 3)
# best_n_cuts <- num_res$results$num_cuts[which.min(num_res$results$BIC)]

# Step 2: Find the exact location of the two best cut-points
cut_res <- find_cutpoint(
  data = lung,
  predictor = "age",
  outcome_time = "time",
  outcome_event = "status",
  method = "genetic",
  num_cuts = 2,
  nmin = 0.1
)

# Step 3: Validate the stability of the cut-points
val_res <- validate_cutpoint(cut_res, num_replicates = 100) # Use more replicates for real analysis

# Plot the final Kaplan-Meier curve and the validation results
plot(cut_res)
plot(val_res)
```

## Getting Started: Full Tutorials

For detailed, step-by-step walkthroughs of these workflows using real-world data, please see the package vignettes, which you can access after installation:

```r
browseVignettes(package = "OptCutR")
```

The vignettes include:

* **Example 1: Diagnostic Cut-point Analysis:** A complete guide to finding and validating a single cut-point for a binary outcome.
* **Example 2: (Placeholder)** *A planned example, e.g., focusing on ecological data.*
* **Example 3: Rapeseed Germination Analysis:** A full three-step workflow for finding and validating multiple cut-points in a survival analysis.

## How to Cite

If you use `OptCutR` in your research, please cite our manuscript:

> Currently in preparation

## Contributing

We welcome contributions and feedback! If you find a bug, have a suggestion, or want to contribute to the code, please open an issue on the [GitHub Issues page](https://github.com/paytonyau/OptCutR/issues).

## License

This package is licensed under the GPL-3 license.
