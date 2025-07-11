# OptCutR: An R Package for Optimal Cut-point Discovery and Validation

<!-- Badges will go here once you set up GitHub Actions and submit to CRAN -->
<!-- e.g., [![CRAN status](https://www.r-pkg.org/badges/version/OptCutR)](https://cran.r-project.org/package=OptCutR) -->
<!-- e.g., [![R-CMD-check](https://github.com/your-username/OptCutR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/your-username/OptCutR/actions/workflows/R-CMD-check.yaml) -->

`OptCutR` provides a comprehensive and flexible workflow to determine, find, and validate optimal cut-points for continuous predictors. It is designed to work with both time-to-event (survival) and binary classification outcomes, making it a versatile tool for researchers in bioinformatics, ecology, medicine, and beyond.

### Key Features

* **Complete Workflow:** Guides you from determining the *best number* of cuts, to finding their *location*, and finally to *validating* the result.
* **Multiple Cut-point Discovery:** Implements a genetic algorithm to find two or more cut-points simultaneously, allowing for stratification into three or more risk groups.
* **Handles Survival and Binary Outcomes:** A unified interface for both time-to-event analysis (using log-rank tests and Cox models) and binary classification (using ROC-based metrics).
* **Robust Methods:** Includes a majority-vote system based on AUC, Youden's J-index, and other metrics to ensure stable results.

## Package Workflow

The `OptCutR` workflow is built around three main functions:


**find_cutpoint_number()**: The first step in an analysis. Use this to get statistical evidence for the most plausible number of cut-points (e.g., 0, 1, or 2+) for your data.

**find_cutpoint()**: The main workhorse function. After you know how many cuts to look for, use this to find their exact locations using either a `"systematic"` search or a `"genetic"` algorithm.

**validate_cutpoint()**: The final step. Use this to assess the stability of your discovered cut-point(s) by running a bootstrap analysis and generating confidence intervals.

The table below outlines the recommended functions and settings for each scenario. 

| Step | Function | Binary Classification | Time-to-Event Analysis (1 Cut-point) | Time-to-Event Analysis (2+ Cut-points) |
| :--- | :--- | :--- | :--- | :--- |
| **1 (Optional)** | `find_cutpoint_number()`  | **Purpose:** Justify dichotomization.<br>**Settings:** `method="systematic"`, `criterion="BIC"` | **Purpose:** Justify dichotomization.<br>**Settings:** `method="systematic"`, `criterion="BIC"` | **Purpose:** Justify using multiple cuts.<br>**Settings:** `method="genetic"`, `max_cuts=2` |
| **2** | `find_cutpoint()` | **Purpose:** Find the best single cut-point.<br>**Settings:** `method="systematic"`, majority vote from ROC metrics. | **Purpose:** Find the best single cut-point.<br>**Settings:** `method="systematic"`, majority vote from survival metrics. | **Purpose:** Find the location of the two best cut-points.<br>**Settings:** `method="genetic"`, `num_cuts=2` |
| **3 (Optional)** | `validate_cutpoint()` | **Purpose:** Assess stability.<br>**Settings:** `num_replicates=500` | **Purpose:** Assess stability.<br>**Settings:** `num_replicates=500` | **Purpose:** Assess stability.<br>**Settings:** `num_replicates=500` |
| **Example** | | [miRNA Diagnostics](https://rpubs.com/payton_yau/OptCutR_1) | [Plant Ecology](https://rpubs.com/payton_yau/OptCutR_2) | [CRC Microbiome](https://rpubs.com/payton_yau/OptCutR_3) |

## Installation

You can install the development version of `OptCutR` from GitHub with:

```r
# install.packages("devtools")
devtools::install_github("paytonyau/OptCutR")
```

*(Once the package is on CRAN, you will be able to install it with `install.packages("OptCutR")`.)*

## Getting Started: A Full Tutorial

For a detailed, step-by-step walkthrough of these workflows using real-world data, please see the "Getting Started" vignette, which you can access after installation:

```r
browseVignettes(package = "OptCutR")
```

## How to Cite

If you use `OptCutR` in your research, please cite our manuscript:

> Currently in preparation 

## Contributing

We welcome contributions and feedback! If you find a bug, have a suggestion, or want to contribute to the code, please open an issue on the [GitHub Issues page](https://github.com/paytonyau/OptCutR/issues).

## License

This package is licensed under the GPL-3 license.