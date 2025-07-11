#' Find the Optimal Number of Cut-points
#'
#' @description
#' Compares models with different numbers of cut-points to determine the most
#' plausible number of cuts for a predictor.
#'
#' @param data A data frame containing the variables.
#' @param predictor The name of the predictor variable, as a string.
#' @param outcome_time The name of the time variable, as a string (for survival).
#' @param outcome_event The name of the event variable, as a string (for survival).
#' @param outcome_binary The name of the binary outcome variable, as a string.
#' @param method The algorithm to use: "systematic" or "genetic".
#' @param criterion The information criterion to use: "BIC" (default) or "AIC".
#' @param max_cuts The maximum number of cut-points to test.
#' @param nmin The minimum number of observations in each group.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores.
#' @param ... Additional arguments passed to the genetic algorithm.
#'
#' @return An object of class `find_cutpoint_number_result`.
#' @importFrom foreach %dopar%
#' @export
find_cutpoint_number <- function(data, predictor,
                                 outcome_time = NULL, outcome_event = NULL,
                                 outcome_binary = NULL,
                                 method = "systematic", criterion = "BIC", max_cuts = 2, nmin = 0.1,
                                 use_parallel = FALSE, ...) {

  # --- 1. Input Validation and Data Prep ---
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("BIC", "AIC"))
  if (is.null(predictor)) stop("A 'predictor' variable must be specified as a string.", call. = FALSE)

  if (!is.null(outcome_time) && !is.null(outcome_event)) {
    analysis_type <- "survival"
    required_vars <- c(predictor, outcome_time, outcome_event)
  } else if (!is.null(outcome_binary)) {
    analysis_type <- "logistic"
    required_vars <- c(predictor, outcome_binary)
  } else {
    stop("You must specify outcome variables as strings.", call. = FALSE)
  }

  if (!all(required_vars %in% names(data))) {
    missing_cols <- required_vars[!required_vars %in% names(data)]
    stop(paste("The following specified columns were not found:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  if (analysis_type == "survival") {
    userdata <- data.frame(time = data[[outcome_time]], event = data[[outcome_event]], factor = data[[predictor]])
  } else {
    userdata <- data.frame(outcome = data[[outcome_binary]], factor = data[[predictor]])
  }
  userdata <- na.omit(userdata)
  n <- nrow(userdata)

  # Process nmin argument
  if (nmin < 1 && nmin > 0) {
    nmin_abs <- ceiling(nmin * n)
    cli::cli_alert_info("Interpreting nmin = {nmin} as a proportion. Minimum group size set to {nmin_abs}.")
    nmin <- nmin_abs
  } else if (nmin >= 1) {
    nmin <- as.integer(nmin)
  } else {
    stop("'nmin' must be a positive number.", call. = FALSE)
  }

  cli::cli_alert_info("Finding optimal number of cuts: method = {.strong {method}}")

  # --- 2. Route to the correct method ---
  if (method == "systematic") {
    results <- .systematic_search(userdata, analysis_type, max_cuts, nmin, criterion)
  } else if (method == "genetic") {
    results <- .genetic_search(userdata, analysis_type, max_cuts, nmin, criterion, ...)
  }

  # --- 3. Unified Post-Processing for AIC/BIC results ---
  if (method %in% c("systematic", "genetic")) {
    min_ic <- min(results[[criterion]], na.rm = TRUE)
    delta_col_name <- paste0("Delta_", criterion)
    weight_col_name <- if (criterion == "AIC") "Akaike_Weight" else "BIC_Weight"

    results[[delta_col_name]] <- results[[criterion]] - min_ic
    exp_delta <- exp(-0.5 * results[[delta_col_name]])
    results[[weight_col_name]] <- exp_delta / sum(exp_delta, na.rm = TRUE)
  }

  output <- list(results = results, parameters = list(method = method, criterion = criterion))
  class(output) <- "find_cutpoint_number_result"

  print(output)
  invisible(output)
}


# --- Internal Helper Functions for find_cutpoint_number ---

.systematic_search <- function(userdata, analysis_type, max_cuts, nmin, criterion) {
    n <- nrow(userdata)

    if (analysis_type == "survival") {
      fit0 <- survival::coxph(survival::Surv(time, event) ~ factor, data = userdata)
      ic0 <- .calc_ic(fit0, 1, n, criterion)
    } else {
      fit0 <- stats::glm(outcome ~ factor, data = userdata, family = stats::binomial())
      ic0 <- .calc_ic(fit0, 2, n, criterion)
    }
    results <- data.frame(num_cuts = 0, IC = ic0)

    for (k_cuts in 1:max_cuts) {
        cli::cli_alert_info("Testing for {k_cuts} cut-point(s)...")
        all_ics_for_k <- c()

        if (k_cuts == 1) {
            grid1 <- userdata$factor[nmin:(n - nmin)]
            for(c1 in grid1) {
                factor_status <- factor(ifelse(userdata$factor <= c1, 0, 1))
                if(min(table(factor_status)) < nmin || nlevels(factor_status) < 2) next
                if(analysis_type == "survival") {
                    fit <- survival::coxph(survival::Surv(time, event) ~ factor_status, data = userdata)
                    all_ics_for_k <- c(all_ics_for_k, .calc_ic(fit, k_cuts, n, criterion))
                } else {
                    fit <- stats::glm(outcome ~ factor_status, data = userdata, family = stats::binomial())
                    all_ics_for_k <- c(all_ics_for_k, .calc_ic(fit, k_cuts + 1, n, criterion))
                }
            }
        } else if (k_cuts == 2) {
            grid1 <- userdata$factor[nmin:(n - 2 * nmin)]
            for(c1 in grid1) {
              start_idx_g2 <- which(userdata$factor > c1)[nmin]
              if(is.na(start_idx_g2)) next
              grid2_end_idx <- n - nmin
              if(start_idx_g2 >= grid2_end_idx) next
              grid2 <- userdata$factor[start_idx_g2:grid2_end_idx]
              for(c2 in grid2){
                if(is.na(c2) || c2 <= c1) next
                factor_status <- as.factor(cut(userdata$factor, breaks=c(-Inf, c1, c2, Inf)))
                if(min(table(factor_status)) < nmin || nlevels(factor_status) < 3) next
                if(analysis_type == "survival") {
                  fit <- survival::coxph(survival::Surv(time, event) ~ factor_status, data = userdata)
                  all_ics_for_k <- c(all_ics_for_k, .calc_ic(fit, k_cuts, n, criterion))
                } else {
                  fit <- stats::glm(outcome ~ factor_status, data = userdata, family = stats::binomial())
                  all_ics_for_k <- c(all_ics_for_k, .calc_ic(fit, k_cuts + 1, n, criterion))
                }
              }
            }
        }
        min_ic_for_k <- if(length(all_ics_for_k) > 0) min(all_ics_for_k, na.rm = TRUE) else NA
        results <- rbind(results, data.frame(num_cuts = k_cuts, IC = min_ic_for_k))
    }
    names(results)[2] <- criterion
    return(results)
}

.genetic_search <- function(userdata, analysis_type, max_cuts, nmin, criterion, ...) {
    if (analysis_type != "survival") stop("The 'genetic' method is only available for survival outcomes.", call. = FALSE)
    if (!exists(".maxloglik", mode = "function")) stop("Helper function .maxloglik not found.", call. = FALSE)
    n <- nrow(userdata)

    fit0 <- survival::coxph(survival::Surv(time, event) ~ factor, data = userdata)
    ic0 <- .calc_ic(fit0, 1, n, criterion)
    results <- data.frame(num_cuts = 0, IC = ic0)

    for (k_cuts in 1:max_cuts) {
        cli::cli_alert_info("Running genetic algorithm for {k_cuts} cut-point(s)...")
        gen_fit <- .maxloglik(target = userdata$factor, numcut = k_cuts, time = userdata$time, censor = userdata$event, confound = NULL, nmin = nmin, ...)
        ic_val <- if(is.null(gen_fit)) NA else .calc_ic(survival::coxph(survival::Surv(time, event) ~ cut(factor, c(-Inf, sort(gen_fit$par[1:k_cuts]), Inf)), data = userdata), k_cuts, n, criterion)
        results <- rbind(results, data.frame(num_cuts = k_cuts, IC = ic_val))
    }
    names(results)[2] <- criterion
    return(results)
}


# --- S3 Methods for Printing and Plotting ---
#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")
  cli::cli_text("Method: {.strong {x$parameters$method}}")
  if(x$parameters$method %in% c("systematic", "genetic")) {
    cli::cli_text("Criterion: {.strong {x$parameters$criterion}}")
  }

  print_df <- x$results
  is_num <- sapply(print_df, is.numeric)
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0) {
      print_df[[weight_col]] <- paste0(round(x$results[[weight_col]] * 100, 1), "%")
  }

  print(print_df, row.names = FALSE)

  if (x$parameters$method %in% c("systematic", "genetic")) {
    best_result <- x$results[which.min(x$results[[x$parameters$criterion]]), ]
    cli::cli_alert_success("\nConclusion: The model with {best_result$num_cuts} cut-point(s) is the most plausible based on {x$parameters$criterion}.")
  }

  cli::cli_text("\nHint: Use plot() on the result object to visualize this table.")
}

#' @export
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results

  if (x$parameters$method %in% c("systematic", "genetic")) {
    y_values <- results[[x$parameters$criterion]]
    y_lab <- x$parameters$criterion
    plot_title <- paste("Model Selection by", x$parameters$criterion)
    plot(results$num_cuts, y_values, type = "b", pch = 19, las = 1,
         xlab = "Number of Cut-points", ylab = y_lab, main = plot_title)
    best_point <- results[which.min(y_values), ]
    points(best_point$num_cuts, min(y_values, na.rm=TRUE), col = "red", pch = 19, cex = 1.5)
  }
}
