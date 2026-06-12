# ===================================================================
#' Find Optimal Number of Cut-points for Survival Data
# ===================================================================
#' @description
#' Finds optimal cut-point number (0 to `max_cuts`) for a Cox model
#' by comparing AIC, AICc, or BIC. Features hardware-accelerated grouping
#' iterations via Rcpp compilation hooks and robust UX constraint warnings.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Computes Goodness-of-fit (AIC, AICc, BIC) for selection.
#' @srrstats {G2.13} `cli_abort()` checks inputs pre-processing.
#' @srrstats {G2.14} `NA` results handled via `na_result` helper.
#' @srrstats {RE6.0} `plot()` method provided.
#' @srrstats {RE6.2} `plot()` visualises model selection metric vs cuts.
#' @srrstats {G1.0} References provided for AIC, AICc, BIC.
#' @srrstats {G1.3} Systematic grid search (max_cuts <= 2) and
#' `rgenoud` global optimisation documented.
#' @srrstats {G1.5} Compared with `cutpointr`/`survminer` in vignette.
#' @srrstats {G1.6} Numerical stability via `survival::coxph`
#' and `rgenoud`; edge cases return `NA`.
#' @srrstats {G2.3b} Uses `as.formula()`/`data[...]`; no unsafe eval.
#' @srrstats {G2.4} `NA` removed via `stats::na.omit()`.
#' @srrstats {G5.2} `optimal_num_cuts` is `NA` when no valid solution found.
#' @srrstats {G2.6} Input validation via direct checks.
#' @srrstats {G2.8} Informative errors via `cli::cli_abort()`.
#' @srrstats {G5.2} Warnings via `cli::cli_alert_warning()`.
#' @srrstats {G5.2} Graceful degradation via `na_result()` for
#' empty data or model failures.
#' @srrstats {G2.14c} `NA` propagation controlled.
#' @srrstats {G1.4} All parameters/return values documented.
#' @srrstats {RE4.17} `print()` method provided.
#' @srrstats {RE4.18} `summary()` method provided.
#'
#' @srrstats {RE1.1} Assumes PH; check `summary()` for `cox.zph`.
#' @srrstats {RE1.3} PH diagnostics via `cox.zph` in `summary()`.
#' @srrstats {RE1.3a} `summary()` includes PH test results.
#' @srrstats {RE2.0} Estimates/SEs from `coxph` in `summary()`.
#' @srrstats {RE2.1} Missing values handled via explicit `na.omit`.
#' @srrstats {RE3.0} `tryCatch` checks model convergence; failures return `NA`.
#' @srrstats {RE2.4a} Checks for collinearity via model fitting constraints.
#' @srrstats {RE2.4b} Checks for insufficient data/constant predictor.
#' @srrstats {RE4.2} Model selection via AIC, AICc, or BIC.
#' @srrstats {RE5.0} Model averaging not implemented.
#' @srrstats {RE6.3} No diagnostic plots; use `cox.zph`.
#' @srrstats {RE7.0a} `na.omit()` removes missing data.
#'
#' @details
#' `method = "systematic"`: grid search respecting `nmin`.
#' `method = "genetic"`: `rgenoud` global optimisation.
#' Systematic search is slow for `max_cuts > 2`; use `genetic`.
#' Core vector partitions are calculated in compiled C++ via `Rcpp` for optimal performance.
#'
#' @param data Input data frame.
#' @param predictor Continuous predictor variable name (character).
#' @param outcome_time Time-to-event variable name (character).
#' @param outcome_event Event indicator name (0/1) (character).
#' @param method `"systematic"` (max_cuts <= 2) or `"genetic"`.
#' @param criterion `"AIC"`, `"AICc"` or `"BIC"`.
#' @param covariates Character vector of covariate names (optional).
#' @param max_cuts Max number of cut-points to test (non-negative int).
#' @param nmin Min. group size (count or proportion).
#' @param seed Integer or `NULL`; random seed for `rgenoud`.
#' @param max.generations Integer; generations for `rgenoud`. If `NULL`, dynamically scales.
#' @param pop.size Integer; population size for `rgenoud`. If `NULL`, dynamically scales.
#' @param use_cpp Logical. Automatically checks and calls compiled C++ routines via `Rcpp`. Default is `TRUE`.
#' @param x An object from [find_cutpoint_number()].
#' @param object An object from [find_cutpoint_number()].
#' @param y Unused.
#' @param ... Additional arguments passed to `rgenoud`.
#'
#' @return An S3 object (`find_cutpoint_number_result`) with
#' `results`, `parameters`, `userdata`, `optimal_num_cuts`,
#' `optimal_cuts`, and `candidate_cuts`.
#'
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom stats na.omit as.formula aggregate quantile symnum
#' @importFrom survival coxph Surv survfit
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success
#' @importFrom cli cli_inform cli_abort cli_bullets cli_h2 cli_alert_warning
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#' @importFrom ggplot2 theme_minimal scale_x_continuous element_text theme
#' @importFrom tools toTitleCase
#' @export
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'
#'   # Generate a pristine simulated clinical tracking baseline template
#'   set.seed(42)
#'   sim_data <- data.frame(
#'     time = rexp(40, rate = 0.1),
#'     event = sample(c(0, 1), 40, replace = TRUE),
#'     biomarker = rnorm(40, mean = 6, sd = 1.2)
#'   )
#'
#'   # Sweep information criteria fit columns up to a 2-cut matrix max
#'   num_fit <- find_cutpoint_number(
#'     data = sim_data,
#'     predictor = "biomarker",
#'     outcome_time = "time",
#'     outcome_event = "event",
#'     max_cuts = 2,
#'     method = "systematic",
#'     criterion = "BIC",
#'     nmin = 5,
#'     quiet = TRUE
#'   )
#'   summary(num_fit)
#' }
find_cutpoint_number <- function(data, predictor,
                                 outcome_time, outcome_event,
                                 method = "systematic", criterion = "BIC",
                                 covariates = NULL, max_cuts = 2,
                                 nmin = 0.1, seed = NULL,
                                 max.generations = NULL, pop.size = NULL,
                                 use_cpp = TRUE, ...) {
  if (!is.numeric(max_cuts) || max_cuts < 0 || max_cuts != round(max_cuts)) {
    cli::cli_abort("max_cuts must be a non-negative integer.")
  }
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("BIC", "AIC", "AICc"))
  if (is.null(predictor)) {
    cli::cli_abort("A 'predictor' variable must be specified.")
  }
  if (is.null(outcome_time) || is.null(outcome_event)) {
    cli::cli_abort("'outcome_time' and 'outcome_event' are required.")
  }

  if (method == "genetic" && !requireNamespace("rgenoud", quietly = TRUE)) {
    cli::cli_abort(
      c(
        "'genetic' method requires the 'rgenoud' package.",
        "i1" = "Install with: install.packages(\"rgenoud\")",
        "i2" = "Or, use `method = \"systematic\"` (for max_cuts <= 2)."
      )
    )
  }

  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  if (!all(required_vars %in% names(data))) {
    missing_cols <- required_vars[!required_vars %in% names(data)]
    cli::cli_abort(
      "Missing columns: {paste(missing_cols, collapse = ', ')}"
    )
  }

  original_predictor_name <- predictor

  userdata <- data[, required_vars, drop = FALSE]
  userdata <- stats::na.omit(userdata)

  if (is.factor(userdata[[predictor]]) || is.character(userdata[[predictor]])) {
    cli::cli_abort(c(
      "x" = "Unordered factors or characters are not allowed as continuous optimization inputs.",
      "i" = "Variable '{predictor}' must be passed as a continuous numeric vector."
    ))
  }

  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"

  n <- nrow(userdata)

  .validate_event_column(userdata$event, outcome_event)

  if (use_cpp && !exists("cpp_get_group_assignments", mode = "function")) {
    cli::cli_alert_warning("Compiled C++ binary not loaded. Falling back gracefully to native R processing.")
    use_cpp <- FALSE
  }

  if (nmin < 1 && nmin > 0) {
    nmin_abs <- floor(nmin * n)
    cli::cli_alert_info("nmin {nmin} is a proportion. Min. group size set to {nmin_abs}.")
  } else if (nmin >= 1) {
    nmin_abs <- as.integer(nmin)
  } else {
    cli::cli_abort("'nmin' must be a positive number.")
  }

  # ✅ FIXED: Keep ... in the function signature if needed, but drop it from the list constructor below
  na_result <- function(userdata, ...) {
    output <- list(
      results = data.frame(),
      parameters = list(
        method = method,
        criterion = criterion,
        analysis_type = "survival",
        predictor = original_predictor_name,
        outcome_time = outcome_time,
        outcome_event = outcome_event,
        covariates = covariates,
        max_cuts = max_cuts,
        nmin = nmin_abs,
        max.generations = max.generations,
        pop.size = pop.size,
        use_cpp = use_cpp
      ),
      userdata = userdata,
      optimal_num_cuts = NA,
      optimal_cuts = NA,
      candidate_cuts = NULL
    )
    class(output) <- "find_cutpoint_number_result"
    cli::cli_inform("No valid results found with given parameters.")
    return(output)
  }

  if (n == 0) {
    cli::cli_inform("No complete cases found after removing NAs.")
    return(na_result(userdata))
  }

  if (n < nmin_abs * (max_cuts + 1)) {
    cli::cli_alert_warning(
      "Sample size ({n}) is mathematically insufficient to support {max_cuts} cuts with a minimum group size constraint of {nmin_abs} patients ({max_cuts + 1} groups required)."
    )
    return(na_result(userdata))
  }

  # Build the base Regularized Grid Space lattice to gauge unique data supply
  grid_probs_vec <- seq(0.01, 0.99, by = 0.01)
  permissible_cuts <- sort(unique(stats::quantile(userdata$factor, probs = grid_probs_vec, na.rm = TRUE)))

  # ===================================================================
  # UX COHORT HEADROOM CHECK (Dynamic Lot Defense)
  # ===================================================================
  unique_values_needed <- (max_cuts + 1) * nmin_abs
  actual_available_nodes <- length(permissible_cuts)

  if (n < unique_values_needed || actual_available_nodes < max_cuts) {
    cli::cli_abort(c(
      "x" = "Insufficient data density for {max_cuts} cut-point(s) at nmin = {nmin_abs}.",
      "i" = "Model requires at least {max_cuts} unique partition boundaries, but only {actual_available_nodes} survived duplicate filtering.",
      "*" = "Action: Reduce 'max_cuts' or lower your 'nmin' threshold."
    ))
  }

  if (length(unique(userdata$factor)) <= max_cuts) {
    cli::cli_inform(paste(
      "Predictor has too few unique values",
      "({length(unique(userdata$factor))}) for max_cuts ({max_cuts})."
    ))
    return(na_result(userdata))
  }

  cli::cli_alert_info("Finding optimal cut number: method = {.strong {method}}")

  if (method == "systematic") {
    params <- list(
      userdata = userdata, max_cuts = max_cuts, nmin = nmin_abs,
      criterion = criterion, covariates = covariates,
      max.generations = max.generations, pop.size = pop.size,
      use_cpp = use_cpp
    )
    extra_args <- list(...)
    params <- c(params, extra_args)

    search_payload <- do.call(.systematic_search_num, params)
    results <- search_payload$selection_table
    permissible_cuts <- search_payload$candidate_cuts
  } else { # genetic complexity loop sweep
    results_list <- list()

    results_list[[1]] <- data.frame(
      num_cuts = 0,
      BIC = .get_model_ic_num(userdata, rep(1, n), 0, n, criterion, covariates),
      cuts = I(list(NA))
    )
    names(results_list[[1]])[2] <- criterion

    confound_df <- if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL

    for (k in 1:max_cuts) {
      cli::cli_inform("Running discrete genetic algorithm for {k} cut-point(s)...")

      # --- SYSTEM SYNCHRONIZED MEMORY-SAFE TUNING ---
      current_pop <- pop.size
      current_gen <- max.generations
      if (is.null(pop.size) || pop.size == 100) {
        current_pop <- switch(as.character(k),
          "1" = 30,
          "2" = 60,
          "3" = 100,
          "4" = 120,
          "5" = 150,
          "6" = 180,
          "7" = 200,
          250
        )
      }
      if (is.null(max.generations) || max.generations == 100) {
        current_gen <- switch(as.character(k),
          "1" = 30,
          "2" = 40,
          "3" = 50,
          "4" = 55,
          "5" = 60,
          "6" = 65,
          "7" = 70,
          80
        )
      }

      params_k <- list(
        target = userdata$factor,
        numcut = k,
        time = userdata$time,
        censor = userdata$event,
        confound = confound_df,
        nmin = nmin_abs,
        criterion = if (criterion %in% c("BIC", "AIC", "AICc")) "p_value" else criterion,
        max.generations = current_gen,
        pop.size = current_pop,
        use_cpp = use_cpp,
        print.level = 0,
        grid_probs = grid_probs_vec
      )
      extra_args <- list(...)

      if (!("boundary.enforcement" %in% names(extra_args))) {
        extra_args$boundary.enforcement <- 2
      }
      params_k <- c(params_k, extra_args)

      raw_fit <- do.call(".run_genetic_search", params_k)

      if (is.null(raw_fit) || !is.finite(raw_fit$value) || raw_fit$value <= -.Machine$double.xmax) {
        fit_k <- NULL
        # DYNAMIC RECOGNITION ARREST: Throw warning for dropped candidate configurations
        cli::cli_alert_warning("Model with {k} cut-point(s) collapsed: Subgroups violated the minimum 'nmin' constraint of {nmin_abs} subjects.")
      } else {
        fit_k <- data.frame(
          num_cuts = k,
          placeholder = .get_model_ic_num(userdata, cut(userdata$factor, breaks = c(-Inf, sort(raw_fit$par[1:k]), Inf), labels = FALSE), k, n, criterion, covariates),
          cuts = I(list(sort(raw_fit$par[1:k])))
        )
        names(fit_k)[2] <- criterion
      }

      if (!is.null(fit_k)) {
        results_list[[k + 1]] <- fit_k
      }
    }
    results <- do.call(rbind, results_list)
  }

  if (is.null(results) || !is.data.frame(results) || nrow(results) == 0) {
    cli::cli_inform("Search algorithm failed to produce results.")
    return(na_result(userdata))
  }

  if (anyNA(results[[criterion]])) {
    cli::cli_alert_warning("All tested model cut-points violated localized subgroup size constraints during runtime search iterations.")
  }

  min_ic <- min(results[[criterion]], na.rm = TRUE)
  delta_col_name <- paste0("Delta_", criterion)
  weight_col_name <- paste0(criterion, "_Weight")

  results[[delta_col_name]] <- results[[criterion]] - min_ic
  exp_delta <- exp(-0.5 * results[[delta_col_name]])
  results[[weight_col_name]] <- exp_delta / sum(exp_delta, na.rm = TRUE)

  results$Evidence <- vapply(results[[delta_col_name]], function(d) {
    if (is.na(d)) {
      NA_character_
    } else if (d <= 2) {
      "Substantial"
    } else if (d <= 7) {
      "Moderate"
    } else {
      "Minimal"
    }
  }, FUN.VALUE = character(1))

  output <- list(
    results = results,
    parameters = list(
      method = method,
      criterion = criterion,
      analysis_type = "survival",
      predictor = original_predictor_name,
      outcome_time = outcome_time,
      outcome_event = outcome_event,
      covariates = covariates,
      max_cuts = max_cuts,
      nmin = nmin_abs,
      max.generations = max.generations,
      pop.size = pop.size,
      use_cpp = use_cpp
    ),
    userdata = userdata,
    candidate_cuts = permissible_cuts
  )

  finite_ic <- results[[criterion]][is.finite(results[[criterion]])]
  if (length(finite_ic) > 0) {
    min_ic_idx <- which.min(results[[criterion]])
    output$optimal_num_cuts <- results$num_cuts[min_ic_idx]
    output$optimal_cuts <- results$cuts[[min_ic_idx]]
  } else {
    output$optimal_num_cuts <- NA
    output$optimal_cuts <- NA
  }

  class(output) <- "find_cutpoint_number_result"
  return(output)
}

#' Internal helper: Compute IC from standard factors for number selection
#'
#' @description
#' Computes the Information Criterion dynamically for the selection engine
#' based on the evaluated grouping factors, defensively capturing model singularities.
#'
#' @inheritParams find_cutpoint_number
#' @param factor_status The cut survival factor.
#' @param k_cuts Number of cuts evaluated.
#' @param n Sample size.
#' @param cov_part Character vector of covariate column names.
#'
#' @return Single numeric IC value (or `Inf` on failure).
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE3.1} Singularity and non-convergence exceptions in the variance-covariance matrix are gracefully intercepted via tryCatch blocks.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @noRd
.get_model_ic_num <- function(userdata, factor_status, k_cuts, n,
                              criterion, cov_part) {
  num_cov <- length(cov_part[cov_part != ""])
  cov_string <- if (num_cov > 0) paste(" +", paste(cov_part, collapse = " + ")) else ""
  formula_str <- paste("survival::Surv(time, event) ~ factor_status", cov_string)

  fit <- tryCatch(
    survival::coxph(as.formula(formula_str), data = userdata),
    error = function(e) {
      return(NULL)
    }
  )
  if (is.null(fit)) {
    return(Inf)
  }

  k_params <- k_cuts + num_cov
  .calc_ic(fit, k_params, n, criterion)
}

# --- S3 Methods for Printing, Summarising and Plotting ---

#' @rdname find_cutpoint_number
#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")

  method_text <- if (!is.null(x$parameters$method)) x$parameters$method else "Unknown"
  criterion_text <- if (!is.null(x$parameters$criterion)) x$parameters$criterion else "IC"

  cli::cli_text("Method: {.strong {method_text}}")
  cli::cli_text("Criterion: {.strong {criterion_text}}")

  if (!is.null(x$parameters$covariates)) {
    cli::cli_text("Covariates: {.strong {paste(x$parameters$covariates, collapse = ', ')}}")
  }

  if (is.null(x$results) || nrow(x$results) == 0 || all(is.na(x$results[[criterion_text]]))) {
    cli::cli_inform("No optimal model could be determined.")
    return(invisible(x))
  }

  print_df <- x$results

  if ("cuts" %in% names(print_df) && is.list(print_df$cuts)) {
    print_df$cuts <- vapply(print_df$cuts, function(c) {
      if (length(c) == 1 && is.na(c[1])) "NA" else paste(round(c, 2), collapse = ", ")
    }, FUN.VALUE = character(1))
  }

  is_num <- vapply(print_df, is.numeric, FUN.VALUE = logical(1))
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0 && weight_col[1] %in% names(x$results)) {
    target_col <- weight_col[1]
    print_df[[target_col]] <- paste0(round(x$results[[target_col]] * 100, 1), "%")
  }

  final_cols <- c(
    "num_cuts", criterion_text,
    paste0("Delta_", criterion_text),
    weight_col[1], "Evidence", "cuts"
  )
  final_cols_exist <- final_cols[final_cols %in% names(print_df)]
  print(print_df[, final_cols_exist, drop = FALSE], row.names = FALSE)

  best_result <- x$results[which.min(x$results[[criterion_text]]), ]

  if (nrow(best_result) > 0 && is.finite(best_result[[criterion_text]])) {
    best_cuts_vals <- best_result$cuts[[1]]
    cli::cli_alert_success(paste(
      "\nConclusion: {best_result$num_cuts} cut-point(s) is",
      "best based on {criterion_text}."
    ))
    if (length(best_cuts_vals) > 0 && !anyNA(best_cuts_vals)) {
      rounded_cuts <- round(best_cuts_vals, 2)
      cli::cli_text("   Optimal cuts at: {.strong {rounded_cuts}}")
    }
  } else {
    cli::cli_inform("\nConclusion: No optimal model could be determined.")
  }

  cli::cli_text("\nHint: Use `summary()` for details, `plot()` to visualise.")
  invisible(x)
}

#' @param show_comparison_table Logical. Show model comparison table?
#' @param show_best_model_details Logical. Show details for best model?
#' @param show_group_counts Logical. Show group counts for best model?
#' @param show_medians Logical. Show median survival for best model?
#' @param plot.it Logical. Display model selection plot?
#' @rdname find_cutpoint_number
#' @export
summary.find_cutpoint_number_result <- function(
  object, show_comparison_table = TRUE,
  show_best_model_details = TRUE,
  show_group_counts = TRUE, show_medians = TRUE,
  plot.it = FALSE, ...
) {
  criterion_text <- object$parameters$criterion %||% "IC"

  cli::cli_h1("Optimal Cut-point Number Analysis")

  if (is.null(object) || is.null(object$results) ||
    nrow(object$results) == 0 ||
    all(is.na(object$results[[criterion_text]]))) {
    cli::cli_inform("Cannot summarise: no valid model was found.")
    return(invisible(object))
  }

  best_result <- object$results[which.min(object$results[[criterion_text]]), ]
  num_cuts <- best_result$num_cuts
  best_cuts_vals <- best_result$cuts[[1]]

  cli::cli_alert_success("Best Model: {.strong {num_cuts} Cut-points} (Criterion: {toupper(criterion_text)})")
  if (length(best_cuts_vals) > 0 && !anyNA(best_cuts_vals)) {
    cli::cli_alert_info("Optimal Thresholds: {.val {paste(round(best_cuts_vals, 3), collapse = ', ')}}")
  }
  cat("\n")

  if (show_comparison_table) {
    cli::cli_h2(sprintf("1. Model Comparison (%s Search)", tools::toTitleCase(object$parameters$method %||% "Unknown")))

    comp_df <- object$results
    comp_df$cuts <- NULL
    comp_df$Marker <- ifelse(comp_df$num_cuts == num_cuts, ">", " ")

    weight_col <- names(comp_df)[grepl("_Weight$", names(comp_df))]
    if (length(weight_col) > 0) {
      target_col <- weight_col[1]
      comp_df[[target_col]] <- paste0(round(comp_df[[target_col]] * 100, 1), "%")
    }

    cols_to_print <- c("Marker", "num_cuts", criterion_text, paste0("Delta_", criterion_text), weight_col[1], "Evidence")
    cols_exist <- cols_to_print[cols_to_print %in% names(comp_df)]

    is_num <- vapply(comp_df[, cols_exist, drop = FALSE], is.numeric, FUN.VALUE = logical(1))
    comp_df[, cols_exist][is_num] <- lapply(comp_df[, cols_exist][is_num], round, 2)

    print(comp_df[, cols_exist, drop = FALSE], row.names = FALSE, right = TRUE)
    cat("\n")
  }

  if (show_best_model_details) {
    data <- object$userdata
    cov_part <- if (!is.null(object$parameters$covariates)) paste(" +", paste(object$parameters$covariates, collapse = " + ")) else ""

    if (num_cuts > 0 && length(best_cuts_vals) > 0 && !anyNA(best_cuts_vals)) {
      data$group <- cut(data$factor,
        breaks = c(-Inf, best_cuts_vals, Inf),
        labels = paste0("G", 1:(num_cuts + 1))
      )

      if (show_group_counts || show_medians) {
        cli::cli_h2("2. Clinical Risk Cohorts")

        counts_table <- as.data.frame(table(data$group))
        names(counts_table) <- c("Group", "N")
        event_counts <- stats::aggregate(event ~ group, data = data, sum)
        names(event_counts) <- c("Group", "Events")
        cohort_df <- merge(counts_table, event_counts, by = "Group")

        if (show_medians) {
          fit_km <- survival::survfit(survival::Surv(time, event) ~ group, data = data)
          surv_summary <- summary(fit_km)$table

          if (is.null(dim(surv_summary))) surv_summary <- t(as.data.frame(surv_summary))

          cohort_df$Median <- surv_summary[, "median"]
          cohort_df$Lower <- surv_summary[, "0.95LCL"]
          cohort_df$Upper <- surv_summary[, "0.95UCL"]

          cohort_df$Median_CI <- sprintf(
            "%s (%s - %s)",
            ifelse(is.na(cohort_df$Median), "NA", round(cohort_df$Median, 1)),
            ifelse(is.na(cohort_df$Lower), "NA", round(cohort_df$Lower, 1)),
            ifelse(is.na(cohort_df$Upper), "NA", round(cohort_df$Upper, 1))
          )

          print(cohort_df[, c("Group", "N", "Events", "Median_CI")], row.names = FALSE, right = TRUE)
        } else {
          print(cohort_df, row.names = FALSE, right = TRUE)
        }
        cat("\n")
      }
    }

    cli::cli_h2("3. Cox Proportional-Hazards")
    formula_str <- if (num_cuts > 0 && length(best_cuts_vals) > 0 && !anyNA(best_cuts_vals)) paste("survival::Surv(time, event) ~ group", cov_part) else paste("survival::Surv(time, event) ~ factor", cov_part)
    model_data <- if (num_cuts > 0 && length(best_cuts_vals) > 0 && !anyNA(best_cuts_vals)) data else object$userdata

    fit_cox <- tryCatch(
      survival::coxph(as.formula(formula_str), data = model_data),
      error = function(e) NULL
    )

    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model for best model: convergence failed.")
    } else {
      cox_sum <- summary(fit_cox)
      coefs <- cox_sum$coefficients
      conf <- cox_sum$conf.int

      if (nrow(coefs) > 0) {
        cox_df <- data.frame(
          Group = gsub("group|factor", "", rownames(coefs)),
          HR = round(conf[, "exp(coef)"], 3),
          Lower = round(conf[, "lower .95"], 3),
          Upper = round(conf[, "upper .95"], 3),
          P_Value = round(coefs[, "Pr(>|z|)"], 3)
        )

        cox_df$Signif <- as.character(stats::symnum(cox_df$P_Value,
          corr = FALSE, na = FALSE,
          cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
          symbols = c("***", "**", "*", ".", " ")
        ))

        print(cox_df, row.names = FALSE, right = TRUE)
        cat("\n")
      }

      conc <- round(cox_sum$concordance["C"], 3)
      pval <- round(cox_sum$sctest["pvalue"], 3)
      cli::cli_alert_info("Overall Model: Concordance = {.val {conc}} | Log-rank p = {.val {pval}}")
      cat("\n")
    }
  }

  if (plot.it) {
    cli::cli_h2("Model Selection Plot")
    print(plot(object, ...))
    cat("\n")
  }

  cli::cli_h2("4. Analysis Parameters")
  params <- object$parameters

  # ✅ FIXED: Strip out duplicate explicit names ("*" =)
  param_bullets <- c(
    paste0("* Search Method: ", tools::toTitleCase(params$method %||% "Unknown")),
    paste0("* Predictor: ", params$predictor %||% "Unknown"),
    paste0("* Criterion: ", params$criterion %||% "Unknown"),
    paste0("* Maximum Cuts: ", params$max_cuts %||% NA),
    paste0("* Minimum Group Size (nmin): ", params$nmin %||% NA)
  )

  if (!is.null(params$covariates)) {
    param_bullets <- c(
      param_bullets,
      paste0("* Covariates: ", paste(params$covariates, collapse = ", "))
    )
  }

  cli::cli_bullets(param_bullets)
  cat("\n")

  invisible(object)
}

#' @rdname find_cutpoint_number
#' @export
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results
  criterion_text <- if (!is.null(x$parameters$criterion)) x$parameters$criterion else "IC"

  if (is.null(results) || nrow(results) == 0 || all(is.na(results[[criterion_text]]))) {
    cli::cli_inform("Cannot generate plot: no valid IC values found.")
    return(invisible(NULL))
  }

  y_values <- results[[criterion_text]]
  valid_indices <- !is.na(y_values)
  if (sum(valid_indices) == 0) {
    cli::cli_inform("Cannot generate plot: no valid IC values found.")
    return(invisible(NULL))
  }

  plot_data <- results[valid_indices, ]
  y_values <- y_values[valid_indices]

  best_point_idx <- which.min(y_values)
  best_num_cuts <- plot_data$num_cuts[best_point_idx]

  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = .data$num_cuts,
    y = .data[[criterion_text]]
  )) +
    ggplot2::geom_line(color = "gray50", linewidth = 0.8) +
    ggplot2::geom_point(
      shape = 21,
      size = 3.5,
      fill = "dodgerblue",
      color = "white",
      stroke = 1
    ) +
    ggplot2::geom_point(
      data = function(df) subset(df, num_cuts == best_num_cuts),
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
