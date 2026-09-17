## Submission summary
This is an update of OptSurvCutR from version 0.10.0 to 0.11.0.

## Test environments
* Local: Windows (R 4.6.1)
* win-builder: R-release, R-devel
* R-hub v2: Linux, Windows, macOS

## R CMD check results
0 errors | 0 warnings | 0 notes

## Changes from previous release
* **Bug fix (Stability Metric Denominator):** Corrected `summary.validate_cutpoint_result()` to compute relative confidence interval widths against the predictor's 10th–90th percentile range by retaining `userdata` in `validate_cutpoint()`, resolving an issue where Stability Tier 1 was unreachable. Removed silent fallback substitutions.
* **Bug fix (Single Cut-point Models):** Single-threshold models are now evaluated on relative width alone, eliminating invalid tier assignments based on interval overlap.
* **API Enhancement:** Added a machine-readable `$stability` list to `summary.validate_cutpoint_result()` for programmatic extraction of diagnostic metrics and tiers.
* **Testing:** Added regression, reachability, and end-to-end workflow integration tests verifying state preservation across `find_cutpoint_number()`, `find_cutpoint()`, and `validate_cutpoint()`.
* **Documentation & Vignettes:** Re-ran case study vignettes to reflect corrected stability metrics and refreshed package documentation and `@srrstats` tags.