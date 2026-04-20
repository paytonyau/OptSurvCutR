## OptSurvCutR 0.2.1

This is a resubmission. In this version, I have addressed the comments from the CRAN reviewer (Konstanze Lauseker):
* Removed single quotes around function and package names in the DESCRIPTION text.
* Replaced instances of `print()` and `cat()` used for informational messages with `message()` in the main functions.
* Removed forced `print()` statements from the end of the main S3 calculation functions, replacing them with standard `return()` calls.
* Bundled a compressed package logo (`man/figures/logo.png`).

## Test environments
* local Windows 11 install, R 4.5.3
* win-builder (devel and release)

## R CMD check results
0 errors | 0 warnings | 0 note

* Note: "New submission" - This is expected as the package is not yet on CRAN.