#' Get a Tidy Summary of Survival Probabilities at Specific Time Points
#'
#' A user-friendly wrapper around `summary.survfit` that automatically handles
#' time unit conversions to provide survival estimates at desired time points
#' (e.g., 1, 3, and 5 years), regardless of the original data's time unit.
#'
#' @param fit An object of class `survfit`, typically created by `survival::survfit()`.
#' @param times A numeric vector of time points at which to get the survival summary.
#' @param data_time_unit The unit of the 'time' variable in the original dataset
#'   used to create the `fit` object. Must be one of "days", "months", or "years".
#' @param summary_time_unit The unit of the `times` argument. Defaults to "years".
#'   Must be one of "days", "months", or "years".
#'
#' @return A tidy data frame summarizing the survival probability, confidence
#'   intervals, and number of subjects at risk for each stratum at each
#'   specified time point.
#'
#' @export
#' @examples
#' \dontrun{
#' library(survival)
#' library(OptCutR)
#'
#' # Load example data where time is in days
#' data(lung)
#' lung$status <- ifelse(lung$status == 2, 1, 0)
#'
#' # Create a survfit object
#' km_fit <- survfit(Surv(time, status) ~ sex, data = lung)
#'
#' # Get the 1-year and 2-year survival probabilities
#' # The function will automatically convert years to days for the summary.
#' surv_summary_at(
#'   fit = km_fit,
#'   times = c(1, 2),
#'   data_time_unit = "days",
#'   summary_time_unit = "years"
#' )
#' }
surv_summary_at <- function(fit, times, data_time_unit, summary_time_unit = "years") {

  # --- 1. Input Validation ---
  if (!inherits(fit, "survfit")) {
    stop("The 'fit' object must be of class 'survfit'.", call. = FALSE)
  }
  valid_units <- c("days", "months", "years")
  if (!data_time_unit %in% valid_units || !summary_time_unit %in% valid_units) {
    stop("'data_time_unit' and 'summary_time_unit' must be one of 'days', 'months', or 'years'.", call. = FALSE)
  }

  # --- 2. Time Unit Conversion ---
  # Define conversion factors relative to 1 day
  conversion_factors <- c(days = 1, months = 30.4375, years = 365.25)

  # Calculate the multiplier to convert user's desired times to the data's scale
  multiplier <- conversion_factors[summary_time_unit] / conversion_factors[data_time_unit]

  # Apply the conversion
  converted_times <- times * multiplier

  # --- 3. Get the Survival Summary ---
  # Use the base summary function with the correctly scaled times
  summary_obj <- summary(fit, times = converted_times, extend = TRUE)

  # --- 4. Format the Output into a Tidy Data Frame ---
  # Check if the model has strata (groups)
  has_strata <- !is.null(summary_obj$strata)

  if (has_strata) {
    # If there are groups, create a data frame for each one and combine
    summary_df <- do.call(rbind, lapply(1:length(summary_obj$strata), function(i) {
      data.frame(
        strata = names(summary_obj$strata)[i],
        time = summary_obj$time[summary_obj$strata == names(summary_obj$strata)[i]],
        n.risk = summary_obj$n.risk[summary_obj$strata == names(summary_obj$strata)[i]],
        n.event = summary_obj$n.event[summary_obj$strata == names(summary_obj$strata)[i]],
        survival = summary_obj$surv[summary_obj$strata == names(summary_obj$strata)[i]],
        lower.ci = summary_obj$lower[summary_obj$strata == names(summary_obj$strata)[i]],
        upper.ci = summary_obj$upper[summary_obj$strata == names(summary_obj$strata)[i]]
      )
    }))
    # Clean up the strata names
    summary_df$strata <- gsub(".*=", "", summary_df$strata)
  } else {
    # If there are no groups, create a simpler data frame
    summary_df <- data.frame(
      time = summary_obj$time,
      n.risk = summary_obj$n.risk,
      n.event = summary_obj$n.event,
      survival = summary_obj$surv,
      lower.ci = summary_obj$lower,
      upper.ci = summary_obj$upper
    )
  }

  # Add the original requested times for clarity
  summary_df$summary_time <- rep(times, length.out = nrow(summary_df))
  summary_df$summary_unit <- summary_time_unit

  # Reorder columns for a clean output
  if (has_strata) {
    summary_df <- summary_df[, c("strata", "summary_time", "summary_unit", "survival", "lower.ci", "upper.ci", "n.risk")]
  } else {
    summary_df <- summary_df[, c("summary_time", "summary_unit", "survival", "lower.ci", "upper.ci", "n.risk")]
  }

  return(summary_df)
}
