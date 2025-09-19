# OptSurvCutR: An R Package for Optimal Cut-point Discovery in Survival Analysis

`OptSurvCutR` provides a comprehensive and flexible workflow to determine, find, and validate optimal cut-points for continuous predictors, with a dedicated focus on **time-to-event (survival)** data. The package's primary goal is to move beyond arbitrary median splits by providing a statistically robust framework to stratify subjects into distinct prognostic groups, with a unique strength in identifying multiple cut-points for complex, non-linear relationships.

### Why OptSurvCutR?

* **Go Beyond the Median Split:** While many tools can find a single cut-point, `OptSurvCutR`'s key strength is its ability to determine the statistically optimal **number** of cut-points, allowing you to uncover complex, non-linear relationships (e.g., U-shaped effects) that a simple dichotomisation would miss.
* **A Complete, End-to-End Workflow:** The package provides a cohesive set of tools that guide you through the entire analytical process: from determining the number of cuts, to finding their precise locations, and finally to validating the stability of your result.
* **Flexible and Powerful Algorithms:** It includes both a systematic grid search for simple cases and a powerful genetic algorithm for efficiently finding multiple cut-points simultaneously.

## Installation

You can install the development version of `OptSurvCutR` from GitHub. Note that the genetic algorithm depends on the `rgenoud` package, which you may need to install first.

```r
# install.packages("remotes")
# install.packages("rgenoud") # Dependency for the genetic algorithm
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Finding Customer Churn Tipping Points

Here is a short example demonstrating the core workflow of `OptSurvCutR` to find optimal `MonthlyCharges` thresholds that predict customer churn.

```r
# Load necessary packages
library(OptSurvCutR)
library(survival)
library(dplyr)

# For this example, we will use a simulated dataset.
set.seed(42)
telco_survival <- data.frame(
  tenure = sample(1:72, 500, replace = TRUE),
  churn_event = sample(0:1, 500, replace = TRUE, prob = c(0.7, 0.3)),
  MonthlyCharges = rnorm(500, 65, 25)
)

# --- 1. Find the optimal NUMBER of cut-points ---
# The BIC suggests 2 cut-points are optimal for this data.
number_result <- find_cutpoint_number(
  data = telco_survival,
  predictor = "MonthlyCharges",
  outcome_time = "tenure",
  outcome_event = "churn_event",
  method = "genetic",
  max_cuts = 3,
  seed = 42,
  quiet = TRUE
)
print(number_result)

# --- 2. Find the optimal VALUE of those cut-points ---
multi_cut_result <- find_cutpoint(
  data = telco_survival,
  predictor = "MonthlyCharges",
  outcome_time = "tenure",
  outcome_event = "churn_event",
  num_cuts = 2, # Use the result from the step above
  method = "genetic",
  seed = 123,
  quiet = TRUE
)
summary(multi_cut_result)

# --- 3. Visualise the Result ---
# The plot reveals three distinct risk groups based on monthly charges.
plot(multi_cut_result, type = "outcome")
```

## Package Workflow

The `OptSurvCutR` workflow is built around three core functions designed to be used in sequence. 

1.  **`find_cutpoint_number()`**: The recommended first step. Use this to get statistical evidence (using BIC, AIC or AICc) for the most plausible number of cut-points for your data.

2.  **`find_cutpoint()`**: The main workhorse function. After you know how many cuts to look for (from Step 1, or if you have a prior hypothesis), use this to find their exact locations.

3.  **`validate_cutpoint()`**: The final step. Use this to assess the stability of your discovered cut-point(s) by running a bootstrap analysis and generating 95% confidence intervals.

## Contributing

We welcome contributions and feedback! Please see the `CONTRIBUTING.md` file for details on how to get involved. If you find a bug, have a suggestion, or want to contribute to the code, please open an issue on the [GitHub Issues page](https://github.com/paytonyau/OptSurvCutR/issues).

Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

## How to Cite

If you use `OptSurvCutR` in your research, please cite our manuscript:

> Currently in preparation

## License

This package is licensed under the GPL-3 license. See the `LICENSE` file for details.
