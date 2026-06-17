#' @keywords internal
#' @srrstats {G1.2} Life Cycle Statement is included in the CONTRIBUTING.md file.
#' @srrstats {G3.0} Verification that statistical algorithms yield expected
#'   outputs has been confirmed via extensive unit tests comparing the internal
#'   C++ routines (`Rcpp`) against native transparent R loop implementations.
#' @srrstats {G5.0} Asserts code pathways verify against defensive data boundaries.
#' @srrstats {G5.1} Software has been tested with varying sample sizes, missing
#'   predictor values (NA handling), and low data-density profiles to ensure
#'   robust handling of boundary conditions.
#' @srrstats {G5.4} Runs low-dimensional matrices to save continuous integration memory.
#' @srrstats {G5.4a} Restricts iterations to single-generation optimisation loops during test sequences.
#' @srrstats {G5.4b} Forces single-core cluster routing as active package defaults during checks.
#' @srrstats {G5.4c} Validates return object class types rather than extracting high-dimensional matrices.
#' @srrstats {G5.6} Fixes seeds to confirm exact numeric search parameter recovery.
#' @srrstats {G5.9} Verifies small input noise values do not shift optimal cut boundaries.
#' @srrstats {RE7.1} For the survival regression workflows, clear descriptions
#'   and diagnostic tools (such as Schoenfeld residual plotting paths) are
#'   provided to explicitly assess the proportional hazards assumption.
"_PACKAGE"

## usethis namespace: start
## usethis namespace: end
NULL

#' NA_standards
#'
#' @description Declarations of non-applicable standards.
#'
#' @srrstatsNA {G2.4e} The package does not implement explicit conversion *from*
#'   factors to other types; it handles factors as grouping variables natively.
#' @srrstatsNA {G3.1a} Arbitrary covariance methods are not supported as the
#'   package relies on standard Cox proportional hazards models (`coxph`) which
#'   handle variance internally.
#' @srrstatsNA {G5.12} Extended tests requiring specific environments are not
#'   implemented; standard tests cover all functionality within CRAN limits.
#' @srrstatsNA {G2.11} This package does not use or support columns with complex
#'   attributes (e.g., 'units') in the predictor or outcome.
#' @srrstatsNA {G2.12} This package does not use or support list-columns in
#'   input data.
#' @srrstatsNA {G2.14b} Missing data is not ignored; it is explicitly removed
#'   (na.omit) to ensure valid survival models.
#' @srrstatsNA {G2.16} NaN/Inf values are handled by the general missing data
#'   removal (na.omit) or trigger errors in model fitting; no specific option to
#'   'keep' them is relevant for Cox models.
#' @srrstatsNA {G3.1} This package does not compute covariance matrices
#'   directly; it relies on 'survival::coxph' for variance estimation.
#' @srrstatsNA {G4.0} This package does not write output to local files.
#' @srrstatsNA {G5.6a} Parameter recovery is tested via statistical significance
#'   and genetic algorithm convergence, not exact parameter value recovery
#'   (which is stochastic for GA).
#' @srrstatsNA {G5.6b} Parameter recovery tests are implicitly covered by the
#'   stochastic genetic algorithm tests.
#' @srrstatsNA {G5.7} Algorithm performance scaling is not explicitly tested as
#'   dataset sizes in this domain are typically small-to-medium (clinical data).
#' @srrstatsNA {G5.8a} Zero-length data is handled by early returns/errors,
#'   validated in G5.2 tests.
#' @srrstatsNA {G5.8b} Unsupported types are handled by input validation (G2.1),
#'   validated in G5.2 tests.
#' @srrstatsNA {G5.8c} All-NA or identical fields are handled by model fitting
#'   checks (coxph failures), validated in G5.2 tests.
#' @srrstatsNA {G5.9a} Trivial noise tests are not implemented as cut-points are
#'   discrete thresholds.
#' @srrstatsNA {G5.9b} Seed tests are implemented (G5.8) to ensure
#'   reproducibility, covering this requirement.
#' @srrstatsNA {RE1.2} Input types are documented in param descriptions;
#'   specific 'numeric' or 'factor' requirements are enforced by validation.
#' @srrstatsNA {RE2.3} Centering or offsetting data is not relevant for finding
#'   cut-points (which are scale-invariant or scale-dependent only on the
#'   predictor's original scale).
#' @srrstatsNA {RE3.1} Convergence warnings from 'coxph' or 'rgenoud' are
#'   captured or suppressed via 'quiet' argument, but fine-grained control is
#'   not exposed.
#' @srrstatsNA {RE3.2} Convergence thresholds are determined by the underlying
#'   'survival' and 'rgenoud' packages; defaults are used.
#' @srrstatsNA {RE3.3} User control of convergence thresholds for internal Cox
#'   models is not exposed; standard defaults are sufficient for cut-point
#'   search.
#' @srrstatsNA {RE4.1} Generating a model object without fitting is not
#'   supported; the package's purpose is to find the fit.
#' @srrstatsNA {RE4.3} Confidence intervals on coefficients are not applicable
#'   as this package returns cut-points, not regression coefficients.
#' @srrstatsNA {RE4.4} Model formula extraction is not relevant for the
#'   cut-point optimisation object.
#' @srrstatsNA {RE4.5} Number of observations is displayed in summary but no
#'   specific `nobs()` S3 method is provided.
#' @srrstatsNA {RE4.6} Variance-covariance matrix (`vcov`) is not applicable to
#'   cut-point optimisation results.
#' @srrstatsNA {RE4.7} Convergence statistics are handled internally by rgenoud
#'   or coxph and not exposed as standard model outputs.
#' @srrstatsNA {RE4.8} Response variables are not stored or accessible via
#'   standard accessor methods on the return object.
#' @srrstatsNA {RE4.9} Modelled values (fitted values) are not returned; the
#'   focus is on group stratification.
#' @srrstatsNA {RE4.10} Residuals are not computed for the cut-point
#'   optimisation itself (though visualised in diagnostic plots).
#' @srrstatsNA {RE4.11} Goodness-of-fit is used for selection (AIC/BIC) but
#'   standard goodness-of-fit accessors are not implemented.
#' @srrstatsNA {RE4.12} Transformations are not part of the model specification
#'   in this context.
#' @srrstatsNA {RE4.13} Predictor variable metadata is not stored in a way
#'   requiring accessors.
#' @srrstatsNA {RE4.14} Forecasting and Extrapolation are out of scope; this
#'   package is for retrospective analysis of existing data.
#' @srrstatsNA {RE4.15} Forecasting error tests are not applicable (no
#'   forecasting capability).
#' @srrstatsNA {RE4.16} Prediction on new groups is not implemented.
#' @srrstatsNA {RE4.18} A summary method is provided, but it does not implement
#'   complex bootstrap summaries for coefficients (bootstrapping is handled by a
#'   separate validation function).
#' @srrstatsNA {RE5.0} Model averaging is not implemented; the package selects a
#'   single optimal solution.
#' @srrstatsNA {RE7.1a} Fitting speed on noiseless data is not a relevant metric
#'   for this algorithm.
#' @srrstatsNA {RE7.2} Return objects retain input metadata (row names not
#'   relevant for aggregate survival data).
#' @srrstatsNA {RE7.3} Accessor methods (coef, etc.) are not implemented (see
#'   RE4.x NAs).
#' @srrstatsNA {RE7.4} Forecast error tests are not applicable.
#' @srrstatsNA {G5.10} Extended tests are not implemented as standard tests
#'   cover the necessary scope.
#' @srrstatsNA {G5.11} Extended tests are not implemented.
#' @srrstatsNA {G5.11a} Extended tests are not implemented.
#'
#' @return NULL. This object is for documentation purposes only.
#' @examples
#' # This is a documentation object for srr standards.
#' NULL
#' @keywords internal
#' @name NA_standards
NULL

# ===================================================================
# Startup Message
# ===================================================================
.onAttach <- function(libname, pkgname) {
  msg <- paste0(
    "== OptSurvCutR v", utils::packageVersion(pkgname), " ==\n",
    "  Docs: <https://github.com/paytonyau/OptSurvCutR>\n",
    "  Paper: Yau, Payton (2025) bioRxiv 10.1101/2025.10.08.681246\n",
    "  Cite: `citation('OptSurvCutR')`\n"
  )
  packageStartupMessage(msg)
}
