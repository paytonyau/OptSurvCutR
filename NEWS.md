# OptSurvCutR 0.1.9.2 (2025-11-25)

## ROPENSCI REVIEW UPDATES (ITERATION 2)

### BREAKING CHANGES
* **Renamed Functions:** `plot_optimization_curve()` $\to$ `plot_optimisation_curve()` (British English) and `plot_schoenfeld()` $\to$ `plot_cutpoint_residuals()` (Namespace conflict).
* **Renamed Argument:** `maxiter` is now `max.generations` in genetic algorithm functions to match `rgenoud`.

### IMPROVEMENTS
* **Explicit Parameters:** Exposed `pop.size` and `max.generations` as explicit arguments in `find_cutpoint()` and `find_cutpoint_number()`.
* **Consistent Logic:** All functions now use `floor()` for `nmin` calculation to ensure identical sample size handling.
* **Better Validation:** `validate_cutpoint()` now accepts `nmin` as a proportion and issues an informative message when using the default 90% buffer (Standard G2.10).
* **Internal Cleanup:** Removed export of helper operator `%||%`.

### DOCUMENTATION & STANDARDS
* **British English:** Standardised all functions and documentation (e.g., "optimise", "summarise").
* **Standards Compliance:** Audited and updated all `@srrstats` tags; moved non-applicable standards to `R/OptSurvCutR-package.R`.
* **Testing:** Added tests for new genetic parameters and updated snapshot tests for S3 methods.

---

# OptSurvCutR 0.1.9.1 (2025-11-20)

## POST-ROPENSCI INITIAL SUBMISSION

* Renamed `plot_diagnostics()` → `plot_schoenfeld()` to avoid name clash with another package.
* Stabilised and accepted all snapshot tests for pathological `find_cutpoint_number_result` S3 methods (now use real criterion `"BIC"`).
* Cleaned DESCRIPTION: removed unused Suggests (`coxphf`, `patchwork`, `pROC`, `tibble`, `srr`) and Remotes (`ropensci-review-tools/srr`).
* Test suite now passes with zero failures/warnings on R CMD check.

---

# OptSurvCutR 0.1.9

## REFACTORING & MAINTENANCE

* **Simplified Parallel Interface:** The `use_parallel` argument was removed from `validate_cutpoint()`. Parallel execution is now controlled **exclusively** by the `n_cores` argument (`n_cores = 1` for sequential, `n_cores > 1` for parallel).
* **Consolidated Validation Logic:** A new internal helper `.validate_event_column()` was created in `R/utils-helpers.R` to centralise the logic for checking that an event column is numeric and contains only 0s and 1s. Both `find_cutpoint_number()` and the internal helper `.validate_data_conditions()` (also in `R/utils-helpers.R`) were updated to call this function, eliminating duplicated code.
* **Simplified Parallel Code:** In `validate_cutpoint()`, the brittle `functions_to_export` block was removed. The function now correctly relies on the `.packages = "OptSurvCutR"` argument in the `foreach` loop, making it easier to maintain.

## BUG FIXES

* **Fixed Parallel Reproducibility:** A bug in `validate_cutpoint()` that prevented true reproducibility for parallel runs (`n_cores > 1`) was fixed. The function now checks for and registers the `{doRNG}` package when a `seed` is provided, ensuring results are identical regardless of the number of cores used. The incorrect `set.seed(i)` call inside the `foreach` loop was removed.
* **Improved S3 Method Robustness:** The `summary.find_cutpoint_number_result()` S3 method was updated to robustly handle `NULL` values in the `object$parameters` list (e.g., `method = NULL`). This prevents potential errors if a result object is created manually or improperly and aligns its behaviour with the `print()` method.

## DOCUMENTATION

* **Clarified Log-Rank Test:** The documentation for `find_cutpoint()` was updated to clarify that when `covariates` are provided, the `"logrank"` criterion is automatically generalized to the more appropriate **Cox score test**.

## IMPROVEMENTS
* **Test Coverage:** Increased test coverage to **88.5%** with new unit tests for parallel reproducibility, S3 method edge cases, and internal validation helpers.
- Planned:  A ROpenSci submission is planned for review.
- Planned:  A JOSS submission is planned post-rOpenSci review.

---

# OptSurvCutR 0.1.8

## CRITICAL BUG FIXES

* **Parallel Validation:** Fixed a critical bug where `validate_cutpoint(use_parallel = TRUE)` would fail. This was caused by helper functions not being exported to the parallel workers.
* **Event Column Validation:** Fixed a major bug where non-numeric `outcome_event` data (e.g., `"0:LIVING"`) caused silent failures. The functions now "fail-fast" with an error, validating the column is numeric with only 0s and 1s.
* **Genetic Algorithm Reporting:** Fixed a bug where `method = "genetic"` would report a failure (`-Inf`) as a successful result. The functions now correctly check for the failure signal and return `NA`.

---

## NEW FEATURES

* **G3.1a Compliance — Schoenfeld Residual Diagnostics**: Added `plot_diagnostics()` to generate publication-ready Schoenfeld residual plots, assessing the proportional hazards (PH) assumption.

---

## IMPROVEMENTS

* **Code Quality (Refactoring):** Refactored `find_cutpoint()` to resolve "high cyclomatic complexity" `NOTE`s by moving validation logic to internal helper functions.
* **Test Coverage:** Increased test coverage by adding new unit tests for error conditions, edge cases, and validation logic.
* **rOpenSci/`pkgcheck` Compliance:**
    * Added examples to all exported functions.
    * Added `@srrstats` tags to satisfy multi-directory requirements.
* **`R CMD check` Compliance:** Addressed all `ERROR`s, `WARNING`s, and `NOTE`s:
    * Added `inst/WORDLIST`, `_PACKAGE` documentation, and fixed empty `Rd` sections.
    * Added `codemeta.json` to `.Rbuildignore`.
* **Code Style:**
    * Reformatted code to follow Tidyverse Style Guide principles (e.g., `<-` for assignment, `"` for strings, consistent spacing around operators).
    * Strictly enforced an 80-character maximum line length for all code, comments, and documentation.
    * Made Roxygen documentation, inline comments, and user-facing messages more concise.
    * Replaced `sapply()` with type-safe `vapply()` in S3 methods to prevent potential bugs.

* **Internal Data:** The `crc_virome` dataset's `status` column was corrected to be numeric 0/1 to match documentation and new validation rules.

---

# OptSurvCutR 0.1.7

## NEW FEATURES
- **Covariate Adjustment Added**: Both `find_cutpoint()` (for `method = "systematic"` and `method = "genetic"`) and `find_cutpoint_number()` now support covariate adjustment via the `covariates` argument. This allows finding optimal cut-points and determining the optimal number of groups while accounting for potential confounders, providing a more robust assessment of a biomarker's independent prognostic value.

## IMPROVEMENTS
- **Major Performance Optimisation**: Significantly optimised performance in `find_cutpoint()` for `criterion = "hazard_ratio"` and `criterion = "p_value"` (both systematic and genetic methods). This was achieved by removing computationally expensive `summary()` calls, extracting coefficients and statistics directly from model objects, and using a fast, manual Likelihood Ratio Test for p-value calculation.
- **Simplified Parallelism**: Removed internal parallel processing (`use_parallel` argument) from `find_cutpoint` and `find_cutpoint_number` systematic search. This prevents potential issues with nested parallel calls and relies on standard external parallelisation approaches, such as the explicit parallel loop within `validate_cutpoint`.
- **Test Suite Overhaul**: Revamped the entire `testthat` suite for improved reliability and robustness. Replaced brittle `expect_snapshot()` tests with more stable checks like `expect_output()`, `expect_s3_class()`, and specific value comparisons. Corrected logic for error and warning expectations and improved mocking for dependency checks. Code coverage increased significantly (for example, to ~86%).
- **User-Friendly `rgenoud` Check**: Added clear, informative error messages using `cli` in `find_cutpoint()` and `find_cutpoint_number()` that trigger immediately if `method = "genetic"` is requested but the suggested `rgenoud` package is not installed, guiding the user on how to install it.
- **Dependency Management**: Moved optional dependencies `broom` (for plotting) and `withr` (for testing) from `Imports` to `Suggests` in the `DESCRIPTION` file, making the core package installation lighter. Removed unused `magrittr` import.
- **Simplified Evidence Labels**: Renamed evidence labels in `find_cutpoint_number()` results for brevity (for example, "Substantial support" -> "Substantial").
- **Code Maintainability**: Centralised all `utils::globalVariables` definitions into `globals.R` to resolve R CMD check `NOTEs` and improve clarity. Removed redundant code.

## BUG FIXES
- **quiet = TRUE Message Fix**: Fixed a critical bug where `find_cutpoint()` and `find_cutpoint_number()` would fail silently (no console message) when `quiet = TRUE` was set. Failure messages are now always printed to the console via internal helpers, regardless of the `quiet` setting, improving user feedback on errors.
- **Genetic Algorithm Edge Cases**: Added more robust input validation and handling within the internal `.run_genetic_search()` function to gracefully manage edge cases like insufficient data variability or non-finite predictor ranges, preventing downstream errors and cryptic warnings.
- **Genetic Algorithm Monitor**: Fixed the internal monitoring function used by the genetic algorithm (`rgenoud::genoud`) to correctly respect the `print.level` argument (controlled indirectly via user functions), ensuring progress updates are displayed or suppressed as intended.
- **NAMESPACE Fix (stats::)**: Resolved `NAMESPACE` errors and related test failures by adding explicit `stats::` calls where needed (for example, for `stats::quantile`, `stats::sd`, `stats::pchisq`) and ensuring correct regeneration of the `NAMESPACE` file via `devtools::document()`.
- **foreach NOTE Fix**: Resolved R CMD check `NOTE` regarding "no visible binding for global".

## IMPROVEMENTS
- Planned:  A ROpenSci submission is planned for review.
- Planned:  A JOSS submission is planned post-rOpenSci review.

---

# OptSurvCutR 0.1.6

## NEW FEATURES
- Added a vignette demonstrating the use of `find_cutpoint()` and `validate_cutpoint()` with TCGA virome data (for example, Alphapapillomavirus as a predictor), guiding users through cut-point optimisation and stability assessment for survival analysis [https://github.com/paytonyau/OptSurvCutR/commit/aa41ca3cb3ff7fdff4cbf6cf8d5de5e4494d3500].
- Introduced comprehensive unit tests using `testthat`, covering core functions (`find_cutpoint()`, `find_cutpoint_number()`, `validate_cutpoint()`) and edge cases like missing data or small sample sizes, with code coverage reporting via `covr` to ensure reliability (>80% coverage) [https://github.com/paytonyau/OptSurvCutR/commit/6579072b087b880b448e101d8eef37f8c4fa5550].

## IMPROVEMENTS
- Optimised the genetic algorithm in `find_cutpoint()` by implementing adaptive `pop.size` (for example, 50 for `num_cuts = 1`) and `max.generations` (for example, 75), reducing runtime by 20–50% for survival datasets while maintaining accuracy for optimal cut-point selection [https://github.com/paytonyau/OptSurvCutR/commit/e0cec20c2a72e48b39e1833742e8d7d829621f39].
- Enhanced error messages in `validate_cutpoint()` to provide specific feedback on bootstrap validation failures, such as insufficient sample sizes or non-converging `coxph` models, improving user debugging experience [https://github.com/paytonyau/OptSurvCutR/commit/92ef362cec36ee9ae976f12a61b65db51cc79d94].
- Added a `pkgdown` GitHub Action to automatically build a package website, improving documentation accessibility, and updated README with badges for build status and code coverage to signal package reliability [https://github.com/paytonyau/OptSurvCutR/commit/38712b8aa3cb919556bdf0e9cba6ca27fda10a60].
- Updated `DESCRIPTION` with corrected URLs, dependency versions, and regenerated Rd files for consistent documentation across all functions [https://github.com/paytonyau/OptSurvCutR/commit/d85d4a14a1a10fb7f68537ea64ef4007c3be465a].

## BUG FIXES
- Fixed NA handling in `find_cutpoint()` to robustly process survival datasets with missing predictor values, preventing errors in `coxph` or `survdiff` model fitting for real-world data like TCGA virome datasets [https://github.com/paytonyau/OptSurvCutR/commit/059a288363c9b7b30272301b709378cb58d76d2b].

---

# OptSurvCutR 0.1.5

## IMPROVEMENTS
- Updated core functions (for example, `find_cutpoint()`, `find_cutpoint_number()`) with improved numerical stability and accuracy for survival model fitting, particularly for genetic algorithm convergence in high-dimensional predictors [https://github.com/paytonyau/OptSurvCutR/commit/9338d0479c34069933562cb6b360428aad9dd6fc].

## BUG FIXES
- Fixed bugs in script handling and input validation, improving reliability for edge cases like small datasets or constant predictors in survival analysis [https://github.com/paytonyau/OptSurvCutR/commit/1de3e1048dcaa33fa8cbb6eab0fa8d89ec5d134c].
- Reverted prior bug fixes to prevent potential regressions, ensuring stable behaviour in `validate_cutpoint()` during bootstrap validation runs [https://github.com/paytonyau/OptSurvCutR/commit/d65fa328634a480a27bde6f64f8d615205e12225].

---

# OptSurvCutR 0.1.0

## NEW FEATURES
- Initial release of `OptSurvCutR` for optimising cut-points in survival analysis.
- Added `find_cutpoint()` to identify optimal cut-points for continuous predictors using systematic or genetic algorithms (via `rgenoud`) with log-rank, p-value, or hazard-ratio criteria.
- Added `find_cutpoint_number()` to select the optimal number of cut-points using AIC, AICc, or BIC.
- Added `validate_cutpoint()` for bootstrap-based stability assessment of cut-points.
- Supports survival models via `survival::coxph` and `survival::survdiff`.
- Includes example usage with simulated churn data, adaptable to TCGA virome datasets (for example, Alphapapillomavirus).

## BUG FIXES
- None (initial release).

## IMPROVEMENTS
- None (initial release).

## DEPRECATIONS
- None (initial release).