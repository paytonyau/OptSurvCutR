

Cran comments · MD
## OptSurvCutR 0.11.1
 
This is a resubmission. In this version, I have addressed the issues
reported by the CRAN incoming checks on the previous submission attempt
(2026-09-18):
 
* Corrected an invalid file URI: `CONTRIBUTING.md` was linked from
  `README.md` using a relative path, which does not resolve on the
  rendered CRAN package page; the file itself had also been committed
  under the wrong case (`contributing.md`). Both are now fixed and the
  link resolves correctly from a clean checkout.
* Replaced the `find_cutpoint_number()` example, which previously took
  10.34s on the win-builder Windows check (flagged as excessive runtime)
  and, independently, exercised a data configuration that returned no
  valid cut-point. The new example runs in under 1 second and
  demonstrates a genuine, well-separated threshold.
* Corrected a stale URL redirect in `NEWS.md` (a FOSDEM 2026 schedule
  link had moved to `archive.fosdem.org`).
In addition, this release fixes two correctness issues found while
preparing the above:
 
* `find_cutpoint_number(method = "systematic")` could silently return a
  0-cut-point model regardless of the data whenever `covariates = NULL`,
  due to an indexing error in the internal log-likelihood-ratio
  calculation. Models fitted with one or more covariates were not
  affected.
* Permutation-adjusted p-values were computed in the wrong direction when
  `criterion = "p_value"` was combined with `n_perm > 0`, understating
  significance. `criterion = "logrank"` and `"hazard_ratio"` were not
  affected.
Both are covered by new regression tests. See `NEWS.md` for full details,
including two further plotting corrections (a mislabelled y-axis default
and a diagnostic panel-labelling issue) bundled into this release.
 
## Test environments
* local Windows 11 install, R 4.6.1
* win-builder (devel and release)
* R CMD check --as-cran, with and without CRAN incoming feasibility checks
## R CMD check results
0 errors | 0 warnings | 0 notes
 
* Note: the CRAN incoming feasibility check caught the CONTRIBUTING.md
  URI issue described above on the first local run of this cycle; it is
  fixed in the version submitted here and the check has been re-run clean.