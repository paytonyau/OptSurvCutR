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
#' @param max_cuts The maximum number of cut-points to test for.
#' @param nmin The minimum number of observations in each group.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores.
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
#'   parameters, and the data used.
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom stats na.omit as.formula pchisq logLik aggregate
#' @importFrom survival coxph Surv survfit
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_warn cli_alert_danger cli_progress_bar
#'  cli_progress_update cli_progress_done
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs theme_minimal scale_x_continuous
#' @export
find_cutpoint_number <- function(data, predictor,
                                 outcome_time, outcome_event,
                                 method = "systematic", criterion = "BIC", max_cuts = 2, nmin = 0.1,
                                 use_parallel = FALSE, seed = NULL, maxiter = 100, ...) {

  # --- 1. Input Validation and Data Prep ---
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("BIC", "AIC", "AICc"))
  if (is.null(predictor)) stop("A 'predictor' variable must be specified as a string.", call. = FALSE)
  if (is.null(outcome_time) || is.null(outcome_event)) stop("Both 'outcome_time' and 'outcome_event' must be specified.", call. = FALSE)

  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  required_vars <- c(predictor, outcome_time, outcome_event)
  if (!all(required_vars %in% names(data))) {
    missing_cols <- required_vars[!required_vars %in% names(data)]
    stop(paste("The following specified columns were not found:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  original_predictor_name <- predictor

  userdata <- data.frame(time = data[[outcome_time]], event = data[[outcome_event]], factor = data[[predictor]])
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

  params <- list(userdata = userdata, max_cuts = max_cuts, nmin = nmin, criterion = criterion,
                 maxiter = maxiter, use_parallel = use_parallel, ...)

  results <- if (method == "systematic") {
    do.call(.systematic_search_num, params)
  } else { # genetic
    do.call(.genetic_search_num, params)
  }

  if (is.null(results) || !is.data.frame(results) || nrow(results) == 0) {
    cli::cli_alert_danger("The search algorithm failed to produce any valid results.")
    return(NULL)
  }

  # --- 3. Unified Post-Processing for AIC/BIC/AICc results ---
  min_ic <- min(results[[criterion]], na.rm = TRUE)
  delta_col_name <- paste0("Delta_", criterion)
  weight_col_name <- paste0(criterion, "_Weight")

  results[[delta_col_name]] <- results[[criterion]] - min_ic
  exp_delta <- exp(-0.5 * results[[delta_col_name]])
  results[[weight_col_name]] <- exp_delta / sum(exp_delta, na.rm = TRUE)

  results$Evidence <- sapply(results[[delta_col_name]], function(d) {
    if (is.na(d)) return(NA_character_)
    if (d <= 2) "Substantial support"
    else if (d <= 7) "Considerably less support"
    else "Essentially no support"
  })

  output <- list(
    results = results,
    parameters = list(
      method = method,
      criterion = criterion,
      analysis_type = "survival",
      predictor = original_predictor_name,
      outcome_time = outcome_time,
      outcome_event = outcome_event
    ),
    userdata = userdata
  )
  class(output) <- "find_cutpoint_number_result"

  print(output)
  invisible(output)
}


# --- Internal Helper Functions ---

.systematic_search_num <- function(userdata, max_cuts, nmin, criterion, use_parallel, ...) {
  if (max_cuts > 2) {
    stop("The 'systematic' method is computationally intensive and is only implemented for max_cuts <= 2.", call. = FALSE)
  }

  n <- nrow(userdata)
  userdata <- userdata[order(userdata$factor), ]

  if (use_parallel) {
    if (!requireNamespace("doParallel", quietly = TRUE)) {
      stop("Package 'doParallel' is required for parallel processing.", call. = FALSE)
    }
    cores <- parallel::detectCores()
    cl <- parallel::makeCluster(cores)
    doParallel::registerDoParallel(cl)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    cli::cli_alert_info("Using {cores} cores for parallel systematic search...")
  } else {
    foreach::registerDoSEQ()
  }

  ic0 <- tryCatch({
    fit0 <- survival::coxph(survival::Surv(time, event) ~ factor, data = userdata)
    .calc_ic(fit0, 1, n, criterion)
  }, error = function(e) {
    cli::cli_warn("Could not calculate IC for the base model (0 cuts). It will be excluded. Error: {e$message}")
    return(NA_real_)
  })

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Testing for {k_cuts} cut-point(s)...")
    best_res_for_k <- list(ic = Inf, cuts = NULL)

    if (k_cuts == 1) {
      grid1 <- unique(userdata$factor[nmin:(n - nmin)])

      res_list <- foreach::foreach(c1 = grid1, .combine = 'rbind', .export = c(".get_model_ic_num", ".calc_ic")) %dopar% {
        factor_status <- factor(ifelse(userdata$factor <= c1, 0, 1))
        if(min(table(factor_status)) < nmin || nlevels(factor_status) < 2) return(NULL)

        current_ic <- .get_model_ic_num(userdata, factor_status, k_cuts, n, criterion)
        if (is.finite(current_ic)) data.frame(ic = current_ic, cuts = c1) else NULL
      }

      if (!is.null(res_list) && nrow(res_list) > 0) {
        best_row <- res_list[which.min(res_list$ic), ]
        best_res_for_k <- list(ic = best_row$ic, cuts = best_row$cuts)
      }

    } else if (k_cuts == 2) {
      grid1 <- unique(userdata$factor[nmin:(nrow(userdata) - 2 * nmin)])

      res_list <- foreach::foreach(c1 = grid1, .combine = 'rbind', .export = c(".get_model_ic_num", ".calc_ic")) %dopar% {
        best_inner_res <- list(ic = Inf, c2 = NA)

        start_idx_g2 <- which(userdata$factor > c1)[nmin]
        if(is.na(start_idx_g2)) return(NULL)
        grid2_end_idx <- n - nmin
        if(start_idx_g2 >= grid2_end_idx) return(NULL)

        grid2 <- unique(userdata$factor[start_idx_g2:grid2_end_idx])
        for(c2 in grid2){
          if(is.na(c2) || c2 <= c1) next
          factor_status <- as.factor(cut(userdata$factor, breaks=c(-Inf, c1, c2, Inf)))
          if(min(table(factor_status)) < nmin || nlevels(factor_status) < 3) next

          current_ic <- .get_model_ic_num(userdata, factor_status, k_cuts, n, criterion)
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

      if (!is.null(res_list) && nrow(res_list) > 0) {
        best_row <- res_list[which.min(res_list$ic), ]
        best_res_for_k <- list(ic = best_row$ic, cuts = c(best_row$c1, best_row$c2))
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

.genetic_search_num <- function(userdata, max_cuts, nmin, criterion, use_parallel, maxiter, ...) {
  n <- nrow(userdata)

  ic0 <- tryCatch({
    fit0 <- survival::coxph(survival::Surv(time, event) ~ factor, data = userdata)
    .calc_ic(fit0, k = 1, n = n, criterion = criterion)
  }, error = function(e) {
    cli::cli_warn("Could not calculate IC for the base model (0 cuts). It will be excluded. Error: {e$message}")
    return(NA_real_)
  })

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Running genetic algorithm for {k_cuts} cut-point(s)...")

    ga_result <- .run_genetic_search(
      target = userdata$factor,
      numcut = k_cuts,
      time = userdata$time,
      censor = userdata$event,
      confound = NULL,
      nmin = nmin,
      criterion = "loglik",
      numgen = maxiter,
      ...
    )

    ic_val <- NA
    cuts_val <- list(NULL)

    if (!is.null(ga_result) && is.finite(ga_result$value)) {
      max_logL <- ga_result$value
      ic_val <- .calc_ic(model = list(loglik = c(NA, max_logL)), k = k_cuts, n = n, criterion = criterion)
      cuts_val <- list(sort(ga_result$par[1:k_cuts]))
    }

    new_row <- data.frame(num_cuts = k_cuts, IC = ic_val)
    new_row$cuts <- I(cuts_val)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

.get_model_ic_num <- function(userdata, factor_status, k_cuts, n, criterion) {
  fit <- tryCatch(survival::coxph(survival::Surv(time, event) ~ factor_status, data = userdata), error = function(e) NULL)
  if(is.null(fit)) return(Inf)
  .calc_ic(fit, k_cuts, n, criterion)
}


# --- S3 Methods for Printing, Summarizing and Plotting ---

#' @rdname find_cutpoint_number
#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")
  cli::cli_text("Method: {.strong {x$parameters$method}}")
  cli::cli_text("Criterion: {.strong {x$parameters$criterion}}")

  print_df <- x$results

  print_df$cuts <- sapply(print_df$cuts, function(c) {
    if (is.null(c)) "NA" else paste(round(c, 2), collapse = ", ")
  })

  is_num <- sapply(print_df, is.numeric)
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0) {
    print_df[[weight_col]] <- paste0(round(x$results[[weight_col]] * 100, 1), "%")
  }

  final_cols <- c("num_cuts", x$parameters$criterion, paste0("Delta_", x$parameters$criterion), weight_col, "Evidence", "cuts")
  final_cols_exist <- final_cols[final_cols %in% names(print_df)]
  print(print_df[, final_cols_exist], row.names = FALSE)

  best_result <- x$results[which.min(x$results[[x$parameters$criterion]]), ]

  if(nrow(best_result) > 0 && is.finite(best_result[[x$parameters$criterion]])){
    best_cuts_vals <- best_result$cuts[[1]]
    cli::cli_alert_success("\nConclusion: The model with {best_result$num_cuts} cut-point(s) is the most plausible based on {x$parameters$criterion}.")
    if (!is.null(best_cuts_vals)) {
      cli::cli_text("  \u2514\u2500 Optimal cuts found at: {.strong {paste(round(best_cuts_vals, 2), collapse = ', ')}}")
    }
  } else {
    cli::cli_alert_danger("\nConclusion: No optimal model could be determined.")
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

  if(is.null(object) || is.null(object$results) || !is.data.frame(object$results)) {
    cli::cli_alert_danger("Cannot generate summary from an invalid or empty result object.")
    return(invisible(NULL))
  }

  if (show_comparison_table) {
    print(object)
  }

  best_result <- object$results[which.min(object$results[[object$parameters$criterion]]), ]

  if(nrow(best_result) == 0 || is.na(best_result[[object$parameters$criterion]])) {
    cli::cli_alert_danger("Cannot generate summary because no valid optimal model was found.")
    return(invisible(NULL))
  }

  if (show_best_model_details) {
    cli::cli_h1("Details for Best Model")

    best_cuts_vals <- best_result$cuts[[1]]
    num_cuts <- best_result$num_cuts

    cli::cli_text("The best model found has {.strong {num_cuts}} cut-point(s).")
    if(!is.null(best_cuts_vals)) {
      cli::cli_text("Cut-point values: {.strong {paste(round(best_cuts_vals, 2), collapse = ', ')}}.")
    }

    data <- object$userdata

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
        fit_km <- survival::survfit(survival::Surv(time, event) ~ group, data = data)
        print(fit_km)
      }
    } else {
      data$group <- data$factor
    }

    cli::cli_h2("Final Cox Proportional-Hazards Model")
    formula_str <- if(num_cuts > 0) "survival::Surv(time, event) ~ group" else "survival::Surv(time, event) ~ factor"
    fit_cox <- survival::coxph(as.formula(formula_str), data = data)
    print(summary(fit_cox))
  }

  if (plot.it) {
    cli::cli_h2("Model Selection Plot")
    print(plot(object, ...))
  }

  invisible(object)
}


#' @rdname find_cutpoint_number
#' @export
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results
  y_values <- results[[x$parameters$criterion]]

  valid_indices <- !is.na(y_values)
  if(sum(valid_indices) == 0) {
    cli::cli_alert_warning("Cannot generate plot because no valid Information Criterion values were calculated.")
    return(invisible(NULL))
  }

  plot_data <- results[valid_indices, ]
  y_values <- y_values[valid_indices]

  best_point_idx <- which.min(y_values)
  best_num_cuts <- plot_data$num_cuts[best_point_idx]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data$num_cuts, y = .data[[x$parameters$criterion]])) +
    ggplot2::geom_line(color = "gray50", linewidth = 0.8) +
    ggplot2::geom_point(shape = 21, size = 3.5, fill = "dodgerblue", color = "white", stroke = 1) +
    ggplot2::geom_point(
      data = ~subset(., num_cuts == best_num_cuts),
      color = "#D55E00", size = 4, shape = 19
    ) +
    ggplot2::scale_x_continuous(breaks = plot_data$num_cuts) +
    ggplot2::labs(
      title = paste("Model Selection by", x$parameters$criterion),
      subtitle = "The best model (lowest value) is highlighted in orange.",
      x = "Number of Cut-points",
      y = x$parameters$criterion
    ) +
    ggplot2::theme_minimal(base_size = 14)

  return(p)
}
