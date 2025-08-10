#' Find the Optimal Number of Cut-points for a Continuous Predictor
#'
#' @description
#' Determines the optimal number of cut-points (from 0 to `max_cuts`) for a
#' continuous variable by comparing models using an information criterion (AIC or BIC).
#' It now also returns the values of the optimal cuts and provides a summary method.
#'
#' @param data A data frame containing the variables.
#' @param predictor The name of the predictor variable, as a string.
#' @param outcome_time The name of the time variable, as a string (for survival).
#' @param outcome_event The name of the event variable, as a string (for survival).
#' @param outcome_binary The name of the binary outcome variable, as a string.
#' @param method The algorithm to use: "systematic" or "genetic".
#' @param criterion The information criterion ("AIC" or "BIC") to use for model selection.
#' @param max_cuts The maximum number of cut-points to test for.
#' @param nmin The minimum number of observations in each group.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores.
#' @param ... Additional arguments passed to the genetic algorithm.
#'
#' @return An object containing the results of the cut-point analysis.
#' @importFrom foreach %dopar%
#' @importFrom stats na.omit as.formula pchisq glm binomial
#' @importFrom survival coxph Surv
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success
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

  # Keep original predictor name for summary
  original_predictor_name <- predictor

  if (analysis_type == "survival") {
    userdata <- data.frame(time = data[[outcome_time]], event = data[[outcome_event]], factor = data[[predictor]])
  } else {
    userdata <- data.frame(outcome = data[[outcome_binary]], factor = data[[predictor]])
  }
  userdata <- na.omit(userdata)
  n <- nrow(userdata)

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
  min_ic <- min(results[[criterion]], na.rm = TRUE)
  delta_col_name <- paste0("Delta_", criterion)
  weight_col_name <- paste0(criterion, "_Weight")

  results[[delta_col_name]] <- results[[criterion]] - min_ic
  exp_delta <- exp(-0.5 * results[[delta_col_name]])
  results[[weight_col_name]] <- exp_delta / sum(exp_delta, na.rm = TRUE)

  ## Improvement: Add a "Strength of Evidence" column
  results$Evidence <- sapply(results[[delta_col_name]], function(d) {
    if (is.na(d)) return(NA)
    if (d <= 2) "Substantial support"
    else if (d <= 7) "Considerably less support"
    else "Essentially no support"
  })

  output <- list(
    results = results,
    parameters = list(
      method = method,
      criterion = criterion,
      analysis_type = analysis_type,
      predictor = original_predictor_name,
      outcome_time = outcome_time,
      outcome_event = outcome_event,
      outcome_binary = outcome_binary
    ),
    userdata = userdata
  )
  class(output) <- "find_cutpoint_number_result"

  print(output)
  invisible(output)
}


# --- Internal Helper Functions ---

.systematic_search <- function(userdata, analysis_type, max_cuts, nmin, criterion) {
  ## Improvement: Add validation for systematic method limitations
  if (max_cuts > 2) {
    stop("The 'systematic' method is computationally intensive and is only implemented for max_cuts <= 2.", call. = FALSE)
  }

  n <- nrow(userdata)
  userdata <- userdata[order(userdata$factor), ]

  if (analysis_type == "survival") {
    fit0 <- survival::coxph(survival::Surv(time, event) ~ factor, data = userdata)
    ic0 <- .calc_ic(fit0, 1, n, criterion)
  } else {
    fit0 <- stats::glm(outcome ~ factor, data = userdata, family = stats::binomial())
    ic0 <- .calc_ic(fit0, 2, n, criterion)
  }
  ## Improvement: Store cut values (NULL for 0 cuts)
  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Testing for {k_cuts} cut-point(s)...")
    best_res_for_k <- list(ic = Inf, cuts = NULL)

    if (k_cuts == 1) {
      grid1 <- unique(userdata$factor[nmin:(n - nmin)])
      for(c1 in grid1) {
        factor_status <- factor(ifelse(userdata$factor <= c1, 0, 1))
        if(min(table(factor_status)) < nmin || nlevels(factor_status) < 2) next

        current_ic <- .get_model_ic(userdata, analysis_type, factor_status, k_cuts, n, criterion)
        if (current_ic < best_res_for_k$ic) {
          best_res_for_k <- list(ic = current_ic, cuts = c1)
        }
      }
    } else if (k_cuts == 2) {
      grid1 <- unique(userdata$factor[nmin:(n - 2 * nmin)])
      for(c1 in grid1) {
        start_idx_g2 <- which(userdata$factor > c1)[nmin]
        if(is.na(start_idx_g2)) next
        grid2_end_idx <- n - nmin
        if(start_idx_g2 >= grid2_end_idx) next

        grid2 <- unique(userdata$factor[start_idx_g2:grid2_end_idx])
        for(c2 in grid2){
          if(is.na(c2) || c2 <= c1) next
          factor_status <- as.factor(cut(userdata$factor, breaks=c(-Inf, c1, c2, Inf)))
          if(min(table(factor_status)) < nmin || nlevels(factor_status) < 3) next

          current_ic <- .get_model_ic(userdata, analysis_type, factor_status, k_cuts, n, criterion)
          if (current_ic < best_res_for_k$ic) {
            best_res_for_k <- list(ic = current_ic, cuts = c(c1, c2))
          }
        }
      }
    }
    min_ic_for_k <- if(is.finite(best_res_for_k$ic)) best_res_for_k$ic else NA
    cuts_for_k <- list(best_res_for_k$cuts)

    new_row <- data.frame(num_cuts = k_cuts, IC = min_ic_for_k)
    new_row$cuts <- I(cuts_for_k)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

.genetic_search <- function(userdata, analysis_type, max_cuts, nmin, criterion, ...) {
  if (analysis_type != "survival") stop("The 'genetic' method is only available for survival outcomes.", call. = FALSE)
  if (!exists(".maxloglik", mode = "function")) stop("Helper function .maxloglik not found. Make sure it's loaded.", call. = FALSE)
  n <- nrow(userdata)

  fit0 <- survival::coxph(survival::Surv(time, event) ~ factor, data = userdata)
  ic0 <- .calc_ic(fit0, 1, n, criterion)
  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Running genetic algorithm for {k_cuts} cut-point(s)...")
    gen_fit <- .maxloglik(target = userdata$factor, numcut = k_cuts, time = userdata$time, censor = userdata$event, confound = NULL, nmin = nmin, ...)

    ic_val <- NA
    cuts_val <- list(NULL)
    if(!is.null(gen_fit)) {
      cuts_val <- list(sort(gen_fit$par[1:k_cuts]))
      fit <- survival::coxph(survival::Surv(time, event) ~ cut(factor, c(-Inf, cuts_val[[1]], Inf)), data = userdata)
      ic_val <- .calc_ic(fit, k_cuts, n, criterion)
    }

    new_row <- data.frame(num_cuts = k_cuts, IC = ic_val)
    new_row$cuts <- I(cuts_val)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

.calc_ic <- function(model, k, n, criterion) {
  logL <- model$loglik[2]
  if (criterion == "BIC") {
    return(-2 * logL + k * log(n))
  } else { # AIC
    return(-2 * logL + 2 * k)
  }
}

.get_model_ic <- function(userdata, analysis_type, factor_status, k_cuts, n, criterion) {
  if(analysis_type == "survival") {
    fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ factor_status, data = userdata), error = function(e) NULL)
    if(is.null(fit)) return(Inf)
    .calc_ic(fit, k_cuts, n, criterion)
  } else {
    fit <- tryCatch(stats::glm(outcome ~ factor_status, data = userdata, family = stats::binomial()), error = function(e) NULL)
    if(is.null(fit)) return(Inf)
    .calc_ic(fit, k_cuts + 1, n, criterion)
  }
}


# --- S3 Methods for Printing, Summarizing and Plotting ---

#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")
  cli::cli_text("Method: {.strong {x$parameters$method}}")
  cli::cli_text("Criterion: {.strong {x$parameters$criterion}}")

  print_df <- x$results

  ## Improvement: Format the cuts nicely for printing
  print_df$cuts <- sapply(print_df$cuts, function(c) {
    if (is.null(c)) "NA" else paste(round(c, 2), collapse = ", ")
  })

  is_num <- sapply(print_df, is.numeric)
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0) {
    print_df[[weight_col]] <- paste0(round(x$results[[weight_col]] * 100, 1), "%")
  }

  # Reorder for clarity
  final_cols <- c("num_cuts", x$parameters$criterion, paste0("Delta_", x$parameters$criterion), weight_col, "Evidence", "cuts")
  print(print_df[, final_cols], row.names = FALSE)

  best_result <- x$results[which.min(x$results[[x$parameters$criterion]]), ]
  best_cuts_vals <- best_result$cuts[[1]]

  cli::cli_alert_success("\nConclusion: The model with {best_result$num_cuts} cut-point(s) is the most plausible based on {x$parameters$criterion}.")

  ## Improvement: Display the optimal cut values
  if (!is.null(best_cuts_vals)) {
    cli::cli_text("  \u2514\u2500 Optimal cuts found at: {.strong {paste(round(best_cuts_vals, 2), collapse = ', ')}}")
  }

  cli::cli_text("\nHint: Use `summary()` for full model details and `plot()` to visualize this table.")
  invisible(x)
}


## Improvement: Add a summary() method
#' @rdname find_cutpoint_number
#' @param object An object of class `find_cutpoint_number_result`.
#' @export
summary.find_cutpoint_number_result <- function(object, ...) {

  best_result <- object$results[which.min(object$results[[object$parameters$criterion]]), ]
  best_cuts_vals <- best_result$cuts[[1]]
  num_cuts <- best_result$num_cuts

  cli::cli_h1("Summary of Best Model")
  cli::cli_text("The best model found has {.strong {num_cuts}} cut-point(s) at {.strong {paste(round(best_cuts_vals, 2), collapse = ', ')}}.")

  data <- object$userdata
  if (num_cuts > 0) {
    data$group <- cut(data$factor, breaks = c(-Inf, best_cuts_vals, Inf),
                      labels = paste0("G", 1:(num_cuts + 1)))
  } else {
    # For the 0-cut model, we just use the continuous predictor
    data$group <- data$factor
  }

  if (object$parameters$analysis_type == "survival") {
    formula_str <- if(num_cuts > 0) "survival::Surv(time, event) ~ group" else "survival::Surv(time, event) ~ factor"
    fit <- survival::coxph(as.formula(formula_str), data = data)
    cli::cli_h2("Final Cox Proportional-Hazards Model")
    print(summary(fit))
  } else { # Logistic
    formula_str <- if(num_cuts > 0) "outcome ~ group" else "outcome ~ factor"
    fit <- stats::glm(as.formula(formula_str), data = data, family = stats::binomial())
    cli::cli_h2("Final Logistic Regression Model")
    print(summary(fit))
  }
  invisible(fit)
}


#' @rdname find_cutpoint_number
#' @param x An object of class `find_cutpoint_number_result`.
#' @param y Unused.
#' @export
#' @importFrom graphics points axis
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results
  y_values <- results[[x$parameters$criterion]]
  y_lab <- x$parameters$criterion
  plot_title <- paste("Model Selection by", x$parameters$criterion)

  plot(results$num_cuts, y_values, type = "b", pch = 19, las = 1,
       xlab = "Number of Cut-points", ylab = y_lab, main = plot_title,
       xaxt = "n")
  axis(1, at = results$num_cuts)

  best_point_idx <- which.min(y_values)
  if(length(best_point_idx) > 0) {
    best_point <- results[best_point_idx, ]
    points(best_point$num_cuts, min(y_values, na.rm=TRUE), col = "red", pch = 19, cex = 1.5)
  }
}
