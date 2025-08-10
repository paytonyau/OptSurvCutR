# OptSurvCutR: An R Package for Optimal Cut-point Discovery in Survival Analysis

<!-- Badges will go here once you set up GitHub Actions and submit to CRAN -->
<!-- e.g., [![CRAN status](https://www.r-pkg.org/badges/version/OptSurvCutR)](https://cran.r-project.org/package=OptSurvCutR) -->
<!-- e.g., [![R-CMD-check](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/paytonyau/OptSurvCutR/actions/workflows/R-CMD-check.yaml) -->

`OptSurvCutR` provides a comprehensive and flexible workflow to determine, find, and validate optimal cut-points for continuous predictors, with a dedicated focus on **time-to-event (survival)** data. It is a versatile tool for any researcher working with this type of data, from bioinformatics and clinical research to economics and sociology. The package's primary goal is to help you find meaningful thresholds that stratify subjects into distinct prognostic or risk groups.

### Key Features

* **Multi-Cut-point Optimisation:** Specifically designed to find **more than one** optimal cut-point, making it ideal for non-linear or U-shaped relationships in your data.

* **Complete Survival Workflow:** Guides you from determining the *best number* of cuts, to finding their *location*, and finally to *validating* the result's stability.

* **Multiple Search Methods:** Find cut-points using either a `"systematic"` brute-force search (for 1-2 cuts) or a fast and powerful `"genetic"` algorithm (for any number of cuts).

* **Multiple Optimisation Criteria:** Go beyond the standard log-rank test. You can find the optimal cut-point by maximising the `"logrank"` statistic, maximising the `"hazard_ratio"`, or minimising the model's `"p_value"`.

* **Evidence-Based Model Selection:** Use information criteria (`"AIC"`, `"AICc"`, or `"BIC"`) to find statistical evidence for the most plausible number of cut-points for your data.

* **Robust Validation:** Includes a built-in `"bootstrapping"` function (`validate_cutpoint`) to assess the stability of your results and generate confidence intervals.

## Package Workflow

The `OptSurvCutR` workflow is built around three core functions designed to be used in sequence:

1. **`find_cutpoint_number()`**: The first step in an analysis. Use this to get statistical evidence for the most plausible number of cut-points (e.g., is 1 cut better than 0 or 2?) for your data.

2. **`find_cutpoint()`**: The main workhorse function. After you know how many cuts to look for, use this to find their exact locations.

3. **`validate_cutpoint()`**: The final step. Use this to assess the stability of your discovered cut-point(s) by running a bootstrap analysis.

The table below outlines the recommended workflow.

| Step | Function | Purpose & Recommended Settings | 
| ----- | ----- | ----- | 
| **1** | `find_cutpoint_number()` | **Purpose:** Justify the number of cut-points.<br>**Settings:** `method="systematic"` or `"genetic"`, `criterion="BIC"` or `"AICc"`. | 
| **2** | `find_cutpoint()` | **Purpose:** Find the location of the cut-point(s).<br>**For 1 cut:** `method="systematic"`, `criterion="logrank"` or `"hazard_ratio"`.<br>**For 2+ cuts:** `method="genetic"`, `num_cuts=2`, `criterion="logrank"`, `use_parallel=TRUE`. | 
| **3** | `validate_cutpoint()` | **Purpose:** Assess the stability of the final cut-point(s).<br>**Settings:** `num_replicates=500`, `use_parallel=TRUE`. | 

## Installation

You can install the development version of `OptSurvCutR` from GitHub. Note that `OptSurvCutR`'s genetic algorithm depends on the `rgenoud` package, which you may need to install first.

```r
# install.packages("devtools")
# install.packages("rgenoud") # Dependency for the genetic algorithm
devtools::install_github("paytonyau/OptSurvCutR")
```

## How to Cite

If you use `OptSurvCutR` in your research, please cite our manuscript:

> Currently in preparation

## Contributing

We welcome contributions and feedback! If you find a bug, have a suggestion, or want to contribute to the code, please open an issue on the [GitHub Issues page](https://github.com/paytonyau/OptSurvCutR/issues).

## License

This package is licensed under the GPL-3 license.
