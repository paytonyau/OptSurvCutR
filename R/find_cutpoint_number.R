#' Find the Optimal Number of Cut-points for Survival Data
#'
#' @description
#' Determines the optimal number of cut-points (from 0 to `max_cuts`) for a
#' continuous variable in a survival context by comparing models using an
#' information criterion (AIC, AICc, or BIC).
#'
#' @param data A data frame containing the variables.
#' @param predictor The name of the predictor variable, as a string.
#' @param outcome_time The name of the time variable, as a string.
#' @param outcome_event The name of the event variable, as a string.
#' @param method The algorithm to use: "systematic" or "genetic".
#' @param criterion The information criterion ("AIC", "AICc", or "BIC") to use for model selection.
#' @param covariates A character vector of covariate names to include in the model.
#' @param max_cuts The maximum number of cut-points to test for.
#' @param nmin The minimum number of observations in each group.
#' @param seed An optional integer for setting the random seed to ensure
#'   reproducible results when using the "genetic" method.
#' @param maxiter The number of generations for the genetic algorithm. Default is 100.
#' @param object An object of class \code{find_cutpoint_number_result} (for S3 methods).
#' @param x An object of class \code{find_cutpoint_number_result} (for S3 methods).
#' @param y Unused.
#' @param ... Additional arguments passed to the genetic algorithm.
#'
#' @return An object of class `find_cutpoint_number_result`. This is a list
#'   containing a data frame of the model comparison results, the analysis
#'   parameters, the data used, optimal_num_cuts, and optimal_cuts.
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom stats na.omit as.formula pchisq logLik aggregate
#' @importFrom survival coxph Surv survfit
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_inform cli_abort
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs theme_minimal scale_x_continuous
#' @export
find_cutpoint_number <- function(data, predictor,
                                 outcome_time, outcome_event,
                                 method = "systematic", criterion = "BIC", covariates = NULL,
                                 max_cuts = 2, nmin = 0.1,
                                 seed = NULL, maxiter = 100, ...) {

  # --- Input Validation ---
  if (!is.numeric(max_cuts) || max_cuts < 0 || max_cuts != round(max_cuts)) {
    cli::cli_abort("max_cuts must be a non-negative integer")
  }
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("BIC", "AIC", "AICc"))
  if (is.null(predictor)) cli::cli_abort("A 'predictor' variable must be specified as a string.")
  if (is.null(outcome_time) || is.null(outcome_event)) cli::cli_abort("Both 'outcome_time' and 'outcome_event' must be specified.")

  # *** NEW: Graceful failure if rgenoud is missing ***
  if (method == "genetic" && !requireNamespace("rgenoud", quietly = TRUE)) {
    cli::cli_abort(
      c("The 'genetic' method requires the 'rgenoud' package.",
        "i" = "Please install it by running: install.packages(\"rgenoud\")",
        "i" = "Alternatively, use `method = \"systematic\"` (for num_cuts <= 2).")
    )
  }

  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  if (!all(required_vars %in% names(data))) {
    missing_cols <- required_vars[!required_vars %in% names(data)]
    cli::cli_abort("The following specified columns were not found: {paste(missing_cols, collapse = ', ')}")
  }

  original_predictor_name <- predictor

  # Rename for internal consistency
  userdata <- data[, required_vars, drop = FALSE]
  userdata <- stats::na.omit(userdata)
  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"

  n <- nrow(userdata)

  # --- Create a standard NA result object ---
  na_result <- function(userdata, method, criterion) {
    output <- list(
      results = data.frame(),
      parameters = list(
        method = method,
        criterion = criterion,
        analysis_type = "survival",
        predictor = original_predictor_name,
        outcome_time = outcome_time,
        outcome_event = outcome_event,
        covariates = covariates, # *** NEW ***
        max_cuts = max_cuts,
        nmin = nmin
      ),
      userdata = userdata,
      optimal_num_cuts = NA,
      optimal_cuts = NA
    )
    class(output) <- "find_cutpoint_number_result"
    if (!is.null(userdata)) cli::cli_inform("No valid results could be determined with the given parameters.")
    return(output)
  }

  # Handle cases with no data after NA removal
  if (n == 0) {
    cli::cli_inform("No complete cases found in the data after removing NAs.")
    return(na_result(userdata, method, criterion))
  }

  if (nmin < 1 && nmin > 0) {
    nmin_abs <- ceiling(nmin * n)
    cli::cli_alert_info("Interpreting nmin = {nmin} as a proportion. Minimum group size set to {nmin_abs}.")
    nmin <- nmin_abs
  } else if (nmin >= 1) {
    nmin <- as.integer(nmin)
  } else {
    cli::cli_abort("'nmin' must be a positive number.")
  }

  # Check for insufficient data
  if (n < nmin * (max_cuts + 1)) {
    cli::cli_inform("Not enough data ({n}) for nmin ({nmin}) and max_cuts ({max_cuts}). Returning empty results.")
    return(na_result(userdata, method, criterion))
  }

  # Check for constant predictor
  if (length(unique(userdata$factor)) <= max_cuts) {
    cli::cli_inform("Predictor has too few unique values ({length(unique(userdata$factor))}) for max_cuts ({max_cuts}). Returning empty results.")
    return(na_result(userdata, method, criterion))
  }

  cli::cli_alert_info("Finding optimal number of cuts: method = {.strong {method}}")

  # Covariates to params
  params <- list(userdata = userdata, max_cuts = max_cuts, nmin = nmin, criterion = criterion,
                 covariates = covariates,
                 maxiter = maxiter, ...)

  results <- if (method == "systematic") {
    do.call(.systematic_search_num, params)
  } else { # genetic
    do.call(.genetic_search_num, params)
  }

  if (is.null(results) || !is.data.frame(results) || nrow(results) == 0) {
    cli::cli_inform("The search algorithm failed to produce any valid results.")
    return(na_result(userdata, method, criterion))
  }

  # --- Unified Post-Processing for AIC/BIC/AICc results ---
  min_ic <- min(results[[criterion]], na.rm = TRUE)
  delta_col_name <- paste0("Delta_", criterion)
  weight_col_name <- paste0(criterion, "_Weight")

  results[[delta_col_name]] <- results[[criterion]] - min_ic
  exp_delta <- exp(-0.5 * results[[delta_col_name]])
  results[[weight_col_name]] <- exp_delta / sum(exp_delta, na.rm = TRUE)

  results$Evidence <- sapply(results[[delta_col_name]], function(d) {
    if (is.na(d)) return(NA_character_)
    if (d <= 2) "Substantial"
    else if (d <= 7) "Moderate"
    else "Minimal"
  })

  # --- Compute optimal_num_cuts and optimal_cuts ---
  output <- list(
    results = results,
    parameters = list(
      method = method,
      criterion = criterion,
      analysis_type = "survival",
      predictor = original_predictor_name,
      outcome_time = outcome_time,
      outcome_event = outcome_event,
      covariates = covariates, # *** NEW ***
      max_cuts = max_cuts,
      nmin = nmin
    ),
    userdata = userdata
  )

  # Add optimal_num_cuts and optimal_cuts
  finite_ic <- results[[criterion]][is.finite(results[[criterion]])]
  if (length(finite_ic) > 0) {
    min_ic_idx <- which.min(results[[criterion]])
    output$optimal_num_cuts <- results$num_cuts[min_ic_idx]
    output$optimal_cuts <- results$cuts[[min_ic_idx]]  # Already a vector or NULL
  } else {
    output$optimal_num_cuts <- NA
    output$optimal_cuts <- NA
  }

  class(output) <- "find_cutpoint_number_result"

  print(output)
  invisible(output)
}

# --- Internal Helper Functions ---

# Added covariates argument
.systematic_search_num <- function(userdata, max_cuts, nmin, criterion, covariates, ...) {
  if (max_cuts > 2) {
    cli::cli_abort("The 'systematic' method is computationally intensive and is only implemented for max_cuts <= 2.")
  }

  n <- nrow(userdata)
  userdata <- userdata[order(userdata$factor), ]

  # Define covariate formula part
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

  # Register sequential backend
  foreach::registerDoSEQ()

  # Base model (0 cuts) - check for continuous factor vs covariates
  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)
  ic0 <- tryCatch({
    fit0 <- survival::coxph(as.formula(base_formula_str), data = userdata)
    .calc_ic(fit0, k = 1 + length(covariates), n = n, criterion = criterion) # k = 1 for factor + covariates
  }, error = function(e) {
    cli::cli_inform("Could not calculate IC for the base model (0 cuts): {e$message}")
    return(NA_real_)
  })

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Testing for {k_cuts} cut-point(s)...")
    best_res_for_k <- list(ic = Inf, cuts = NULL)

    if (k_cuts == 1) {
      grid1 <- unique(userdata$factor[nmin:(n - nmin)])
      if (length(grid1) == 0) {
        cli::cli_inform("Not enough data ({n}) for nmin ({nmin}) and {k_cuts} cut(s). Skipping.")
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      res_list <- foreach::foreach(c1 = grid1, .combine = 'rbind', .export = c(".get_model_ic_num", ".calc_ic")) %do% {
        factor_status <- factor(ifelse(userdata$factor <= c1, 0, 1))
        if (min(table(factor_status)) < nmin || nlevels(factor_status) < 2) return(NULL)

        # *** NEW: Pass cov_part ***
        current_ic <- .get_model_ic_num(userdata, factor_status, k_cuts, n, criterion, cov_part)
        if (is.finite(current_ic)) data.frame(ic = current_ic, cuts = c1) else NULL
      }

      if (is.null(res_list) || nrow(res_list) == 0) {
        cli::cli_inform("No valid cut-points found for {k_cuts} cut(s) due to model failures or constraints.")
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      best_row <- res_list[which.min(res_list$ic), ]
      best_res_for_k <- list(ic = best_row$ic, cuts = best_row$cuts)

    } else if (k_cuts == 2) {
      grid1 <- unique(userdata$factor[nmin:(nrow(userdata) - 2 * nmin)])
      if (length(grid1) == 0) {
        cli::cli_inform("Not enough data ({n}) for nmin ({nmin}) and {k_cuts} cut(s). Skipping.")
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      res_list <- foreach::foreach(c1 = grid1, .combine = 'rbind', .export = c(".get_model_ic_num", ".calc_ic")) %do% {
        best_inner_res <- list(ic = Inf, c2 = NA)

        start_idx_g2 <- which(userdata$factor > c1)[nmin]
        if (is.na(start_idx_g2)) return(NULL)
        grid2_end_idx <- n - nmin
        if (start_idx_g2 >= grid2_end_idx) return(NULL)

        grid2 <- unique(userdata$factor[start_idx_g2:grid2_end_idx])
        if (length(grid2) == 0) return(NULL)

        for (c2 in grid2) {
          if (is.na(c2) || c2 <= c1) next
          factor_status <- as.factor(cut(userdata$factor, breaks = c(-Inf, c1, c2, Inf)))
          if (min(table(factor_status)) < nmin || nlevels(factor_status) < 3) next

          # Pass cov_part
          current_ic <- .get_model_ic_num(userdata, factor_status, k_cuts, n, criterion, cov_part)
          if (current_ic < best_inner_res$ic) {
            best_inner_res <- list(ic = current_ic, c2 = c2)
          }
        }

        if (is.finite(best_inner_res$ic)) {
          data.frame(ic = best_inner_res$ic, c1 = c1, c2 = best_inner_res$c2)
        } else {
          NULL
        }
      }

      if (is.null(res_list) || nrow(res_list) == 0) {
        cli::cli_inform("No valid cut-points found for {k_cuts} cut(s) due to model failures or constraints.")
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      best_row <- res_list[which.min(res_list$ic), ]
      best_res_for_k <- list(ic = best_row$ic, cuts = c(best_row$c1, best_row$c2))
    }

    min_ic_for_k <- if (is.finite(best_res_for_k$ic)) best_res_for_k$ic else NA_real_
    cuts_for_k <- list(best_res_for_k$cuts)

    new_row <- data.frame(num_cuts = k_cuts, IC = min_ic_for_k)
    new_row$cuts <- I(cuts_for_k)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

# Added covariates argument
.genetic_search_num <- function(userdata, max_cuts, nmin, criterion, covariates, maxiter, ...) {
  n <- nrow(userdata)

  # Define covariate formula part
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""
  num_cov <- length(covariates)

  # Base model (0 cuts) - check for continuous factor vs covariates
  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)
  ic0 <- tryCatch({
    fit0 <- survival::coxph(as.formula(base_formula_str), data = userdata)
    .calc_ic(fit0, k = 1 + num_cov, n = n, criterion = criterion) # k = 1 for factor + covariates
  }, error = function(e) {
    cli::cli_inform("Could not calculate IC for the base model (0 cuts): {e$message}")
    return(NA_real_)
  })

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Running genetic algorithm for {k_cuts} cut-point(s)...")

    ga_result <- tryCatch({
      .run_genetic_search(
        target = userdata$factor,
        numcut = k_cuts,
        time = userdata$time,
        censor = userdata$event,
        confound = if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL, # *** NEW ***
        nmin = nmin,
        criterion = "loglik", # Always use loglik for IC calculation
        numgen = maxiter,
        ...
      )
    }, error = function(e) {
      cli::cli_inform("Genetic algorithm failed for {k_cuts} cut(s): {e$message}")
      return(NULL)
    })

    ic_val <- NA_real_
    cuts_val <- list(NULL)

    if (!is.null(ga_result) && is.finite(ga_result$value)) {
      max_logL <- ga_result$value
      # k = num_cuts (for k+1 groups, so k params) + num_cov
      k_params <- k_cuts + num_cov
      ic_val <- .calc_ic(model = list(loglik = c(NA, max_logL)), k = k_params, n = n, criterion = criterion)
      cuts_val <- list(sort(ga_result$par[1:k_cuts]))
    } else {
      cli::cli_inform("No valid cut-points found for {k_cuts} cut(s) due to genetic algorithm failure or constraints.")
    }

    new_row <- data.frame(num_cuts = k_cuts, IC = ic_val)
    new_row$cuts <- I(cuts_val)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

# Added cov_part argument
.get_model_ic_num <- function(userdata, factor_status, k_cuts, n, criterion, cov_part) {
  num_cov <- length(cov_part[cov_part != ""]) # Count covariates

  # Build formula with covariates
  formula_str <- paste("survival::Surv(time, event) ~ factor_status", cov_part)

  fit <- tryCatch(
    survival::coxph(as.formula(formula_str), data = userdata),
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(fit)) return(Inf)

  # k = k_cuts (for k+1 groups, so k params) + num_cov
  k_params <- k_cuts + num_cov
  .calc_ic(fit, k_params, n, criterion)
}


# --- S3 Methods for Printing, Summarizing and Plotting ---

#' @rdname find_cutpoint_number
#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")

  method_text <- if (!is.null(x$parameters$method)) x$parameters$method else "Unknown"
  criterion_text <- if (!is.null(x$parameters$criterion)) x$parameters$criterion else "IC"

  cli::cli_text("Method: {.strong {method_text}}")
  cli::cli_text("Criterion: {.strong {criterion_text}}")

  # Print covariates
  if (!is.null(x$parameters$covariates)) {
    cli::cli_text("Covariates: {.strong {paste(x$parameters$covariates, collapse = ', ')}}")
  }

  if (is.null(x$results) || nrow(x$results) == 0 || all(is.na(x$results[[criterion_text]]))) {
    cli::cli_inform("No optimal model could be determined.")
    return(invisible(x))
  }

  print_df <- x$results

  if ("cuts" %in% names(print_df) && is.list(print_df$cuts)) {
    print_df$cuts <- sapply(print_df$cuts, function(c) {
      if (is.null(c)) "NA" else paste(round(c, 2), collapse = ", ")
    })
  }

  is_num <- sapply(print_df, is.numeric)
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0 && weight_col %in% names(x$results)) {
    print_df[[weight_col]] <- paste0(round(x$results[[weight_col]] * 100, 1), "%")
  }

  final_cols <- c("num_cuts", criterion_text, paste0("Delta_", criterion_text), weight_col, "Evidence", "cuts")
  final_cols_exist <- final_cols[final_cols %in% names(print_df)]
  print(print_df[, final_cols_exist, drop = FALSE], row.names = FALSE)

  best_result <- x$results[which.min(x$results[[criterion_text]]), ]

  if (nrow(best_result) > 0 && is.finite(best_result[[criterion_text]])) {
    best_cuts_vals <- best_result$cuts[[1]]
    cli::cli_alert_success("\nConclusion: The model with {best_result$num_cuts} cut-point(s) is the most plausible based on {criterion_text}.")
    if (!is.null(best_cuts_vals)) {
      cli::cli_text("  \u2514\u2500 Optimal cuts found at: {.strong {paste(round(best_cuts_vals, 2), collapse = ', ')}}")
    }
  } else {
    cli::cli_inform("\nConclusion: No optimal model could be determined.")
  }

  cli::cli_text("\nHint: Use `summary()` for full model details and `plot()` to visualize this table.")
  invisible(x)
}

#' @param show_comparison_table Logical. If TRUE, shows the model comparison table.
#' @param show_best_model_details Logical. If TRUE, shows full details for the best model.
#' @param show_group_counts Logical. If TRUE, shows group counts for the best model.
#' @param show_medians Logical. If TRUE, shows median survival for the best model.
#' @param plot.it Logical. If TRUE, displays the model selection plot.
#' @rdname find_cutpoint_number
#' @export
summary.find_cutpoint_number_result <- function(object, show_comparison_table = TRUE, show_best_model_details = TRUE,
                                                show_group_counts = TRUE, show_medians = TRUE, plot.it = FALSE, ...) {

  criterion_text <- if (!is.null(object$parameters$criterion)) object$parameters$criterion else "IC"

  cli::cli_h1("Optimal Cut-point Number Analysis ({tools::toTitleCase(object$parameters$method)})")

  if (is.null(object) || is.null(object$results) || nrow(object$results) == 0 || all(is.na(object$results[[criterion_text]]))) {
    cli::cli_inform("Cannot generate summary because no valid optimal model was found.")
    return(invisible(object))
  }

  if (show_comparison_table) {
    print(object)
  }

  if (show_best_model_details) {
    cli::cli_h1("Details for Best Model")

    best_result <- object$results[which.min(object$results[[criterion_text]]), ]
    best_cuts_vals <- best_result$cuts[[1]]
    num_cuts <- best_result$num_cuts

    cli::cli_text("The best model found has {.strong {num_cuts}} cut-point(s).")
    if (!is.null(best_cuts_vals)) {
      cli::cli_text("Cut-point values: {.strong {paste(round(best_cuts_vals, 2), collapse = ', ')}}.")
    }

    data <- object$userdata

    # *** NEW: Handle covariates in formula ***
    cov_part <- if (!is.null(object$parameters$covariates)) {
      paste(" +", paste(object$parameters$covariates, collapse = " + "))
    } else {
      ""
    }

    if (num_cuts > 0) {
      data$group <- cut(data$factor, breaks = c(-Inf, best_cuts_vals, Inf),
                        labels = paste0("G", 1:(num_cuts + 1)))

      if (show_group_counts) {
        cli::cli_h2("Group Counts")
        counts_table <- as.data.frame(table(data$group))
        names(counts_table) <- c("Group", "N")
        event_counts <- stats::aggregate(event ~ group, data = data, sum)
        names(event_counts) <- c("Group", "Events")
        print(merge(counts_table, event_counts, by = "Group"))
      }

      if (show_medians) {
        cli::cli_h2("Median Survival by Group")
        # Note: Median survival is unadjusted (doesn't use covariates)
        fit_km <- survival::survfit(survival::Surv(time, event) ~ group, data = data)
        print(fit_km)
      }
    }

    cli::cli_h2("Final Cox Proportional-Hazards Model")
    formula_str <- if (num_cuts > 0) {
      paste("survival::Surv(time, event) ~ group", cov_part)
    } else {
      paste("survival::Surv(time, event) ~ factor", cov_part)
    }

    model_data <- if (num_cuts > 0) data else object$userdata

    fit_cox <- tryCatch(survival::coxph(as.formula(formula_str), data = model_data), error = function(e) NULL)
    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model for best model: model convergence failed.")
    } else {
      print(summary(fit_cox))
    }
  }

  if (plot.it) {
    cli::cli_h2("Model Selection Plot")
    print(plot(object, ...))
  }

  cli::cli_h1("Analysis Parameters")
  params <- object$parameters
  param_bullets <- c(
    "*" = "Search Method: {tools::toTitleCase(params$method)}",
    "*" = "Predictor: {params$predictor}",
    "*" = "Criterion: {params$criterion}",
    "*" = "Maximum Cuts: {params$max_cuts}",
    "*" = "Minimum Group Size (nmin): {params$nmin}"
  )
  # *** NEW: Add covariates to summary ***
  if (!is.null(params$covariates)) {
    param_bullets <- c(param_bullets, "*" = "Covariates: {paste(params$covariates, collapse = ', ')}")
  }

  cli::cli_bullets(param_bullets)

  invisible(object)
}

#' @rdname find_cutpoint_number
#' @export
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results
  criterion_text <- if (!is.null(x$parameters$criterion)) x$parameters$criterion else "IC"

  if (is.null(results) || nrow(results) == 0 || all(is.na(results[[criterion_text]]))) {
    cli::cli_inform("Cannot generate plot because no valid Information Criterion values were calculated.")
    return(invisible(NULL))
  }

  y_values <- results[[criterion_text]]
  valid_indices <- !is.na(y_values)
  if (sum(valid_indices) == 0) {
    cli::cli_inform("Cannot generate plot because no valid Information Criterion values were calculated.")
    return(invisible(NULL))
  }

  plot_data <- results[valid_indices, ]
  y_values <- y_values[valid_indices]

  best_point_idx <- which.min(y_values)
  best_num_cuts <- plot_data$num_cuts[best_point_idx]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$num_cuts,
                                               y = .data[[criterion_text]])) +
    ggplot2::geom_line(color = "gray50", linewidth = 0.8) +
    ggplot2::geom_point(shape = 21,
                        size = 3.5,
                        fill = "dodgerblue",
                        color = "white",
                        stroke = 1) +
    ggplot2::geom_point(
      data = ~subset(., num_cuts == best_num_cuts),
      color = "#D55E00", size = 4, shape = 19
    ) +
    ggplot2::scale_x_continuous(breaks = plot_data$num_cuts) +
    ggplot2::labs(
      title = paste("Model Selection by", criterion_text),
      subtitle = "Best model (lowest) in orange",
      x = "Number of Cut-points",
      y = criterion_text
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 10))

  return(p)
}

