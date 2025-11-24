# ===================================================================
#' Global Variable Definitions for NSE
# ===================================================================

#' @description
#' Declares variable names used in NSE contexts (e.g., `ggplot2::aes()`)
#' to suppress `R CMD check` "no visible binding" notes.
#'
#' Not exported; internal use only.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @noRd
NULL

# NSE global variable definitions
utils::globalVariables(c(
  # --- General / Legacy ---
  "density", "Value", "everything", "c1",
  # --- Effect size & waterfall plots ---
  "Cut1", "HR", "HR_low", "HR_up",
  "OR", "OR_low", "OR_up",
  "classified_group", "is_correct", "outcome", "patient_id",
  # --- Cut-point plotting ---
  "cut1", "stat", "factor", "original_cut", "num_cuts",
  # --- Confidence intervals ---
  "ci_lower", "ci_upper"
))
