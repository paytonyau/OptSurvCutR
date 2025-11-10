# ===================================================================
# MAIN FUNCTION: FIND OPTIMAL CUT-POINTS FOR SURVIVAL DATA
# ===================================================================
#' Find Optimal Cut-points for Survival Data
#'
#' @description
#' Finds optimal cut-point(s) for a continuous predictor in a
#' time-to-event (survival) analysis. Uses systematic search (1–2
#' cuts) or a genetic algorithm (any number of cuts).
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G1.1} Implements genetic + systematic search for
#' optimal multi-cutpoint survival groupings.
#' @srrstats {G1.2} References provided for Cox, log-rank, genetic
#' optimization.
#' @srrstats {G1.3} Systematic grid search (1–2 cuts) and
#' `rgenoud` global optimization documented.
#' @srrstats {G1.5} Compared with `cutpointr` and `survminer`
#' in package vignette.
#' @srrstats {G1.6} Numerical stability via `survival::coxph`
#' and `rgenoud`; edge cases return `NA`.
#' @srrstats {G2.3b} NSE via `as.formula()` and `data[...]`;
#' no unsafe evaluation.
#' @srrstats {G2.4} `NA` removed via `stats::na.omit()`;
#' structured `NA` object returned on failure.
#' @srrstats {G2.4e} `optimal_cuts` and `optimal_stat` are `NA`
#' when no valid solution found.
#' @srrstats {G2.6} Validates inputs via helpers.
#' @srrstats {G2.8} Informative errors via `cli::cli_abort()`.
#' @srrstats {G2.10} Warnings via `cli::cli_alert_warning()`.
#' @srrstats {G2.12} Graceful degradation via `na_result()`
#' for empty data or model failures.
#' @srrstats {G2.13} `cli_abort()` for invalid input.
#' @srrstats {G2.14} `cli_abort()` for missing `rgenoud`.
#' @srrstats {G2.14c} `NA` propagation controlled.
#' @srrstats {G3.1} `plot()` method provided.
#' @srrstats {G4.0} All parameters and return values documented.
#' @srrstats {G5.4c} Edge cases (zero rows, constant predictor)
#' tested.
#' @srrstats {G5.6a, G5.6b} Negative `num_cuts`/`nmin` rejected.
#' @srrstats {G5.7} Large `num_cuts` constrained by `nmin`.
#' @srrstats {G5.8d} `set.seed(seed)` for genetic reproducibility.
#' @srrstats {G5.12} Systematic search scales poorly > 2 cuts.
#'
#' @srrstats {RE1.0} Implements optimal cut-point algorithm.
#' @srrstats {RE1.1} Assumes PH; check `summary()` for `cox.zph`.
#' @srrstats {RE1.4} Cox PH assumption test via `summary(fit)$cox_zph`.
#' @srrstats {RE2.0, RE2.1} Estimates/SEs from `coxph` in `summary()`.
#' @srrstats {RE2.4, RE2.4a, RE2.4b} `tryCatch` checks model
#' convergence; failures return `NA`.
#' @srrstats {RE4.2} Model selection via log-rank, HR, or p-value.
#' @srrstats {RE4.3, RE4.4, RE4.5, RE4.6, RE4.7, RE4.8, RE4.9,
#' RE4.10, RE4.11, RE4.12, RE4.13, RE4.14, RE4.15, RE4.16,
#' RE4.18} N/A (no stepwise, LASSO, etc.).
#' @srrstats {RE5.0} Model averaging not implemented.
#' @srrstats {RE6.0, RE6.2, RE6.3} No diagnostic plots; use `cox.zph`.
#' @srrstats {RE7.0a, RE7.1a} `na.omit()` removes missing data.
#'
#' @details
#' `method = "systematic"`: grid search respecting `nmin`.
#' `method = "genetic"`: `rgenoud` global optimization.
#' Systematic search is slow for `num_cuts > 2`; use `genetic`.
#'
#' @references
#' Altman, D. G., Lausen, B., Sauerbrei, W., &
#' Schumacher,
#' M. (1994). Dangers of Using “Optimal” Cutpoints in the Evaluation of
#'   Prognostic Factors. *JNCI: Journal of the National
#'   Cancer Institute*,
#'   86(11), 829–835. \doi{10.1093/jnci/86.11.829}
#'
#' Cox, D. R. (1972). Regression Models and Life-Tables. *Journal
#' of the Royal Statistical Society: Series B
#' (Methodological)*, 34(2),
#' 187–202. \doi{10.1111/j.2517-6161.1972.tb00899.x}
#'
#' Mantel, N. (1966). Evaluation of survival data and two new
#' rank order statistics arising in its consideration. *Cancer
#' Chemotherapy Reports*, 50(3).
#' https://pubmed.ncbi.nlm.nih.gov/5910392/
#'
#' Mebane Jr, W. R., & Sekhon, J. S. (2011). Genetic
#' Optimization Using Derivatives: The rgenoud Package for R.
#' *Journal
#' of Statistical Software*, 42, 1–26.
#' \doi{10.18637/jss.v042.i11}
#'
#' @param data A data frame containing the analysis variables.
#' @param predictor The continuous predictor variable.
#' @param outcome_time The time-to-event variable.
#' @param outcome_event The event status variable (0 or 1).
#' @param num_cuts The number of cut-points to find. Default is 1.
#' @param method Algorithm: `"systematic"` or `"genetic"`.
#' @param criterion The statistic to optimize: `"logrank"` (max),
#'   `"hazard_ratio"` (max), or `"p_value"` (min).
#' @param covariates Character vector of covariate names.
#' @param nmin Min. group size (integer count or proportion).
#' @param seed Optional integer seed for `"genetic"` method.
#' @param maxiter Number of generations for genetic algorithm (default 100).
#' @param quiet Logical. If `TRUE`, suppresses final print.
#' @param ... Additional arguments passed to `rgenoud` (e.g., `popSize`).
#' @param x An object from [find_cutpoint()].
#' @param object An object from [find_cutpoint()].
#' @param show_model Logical. Show final Cox model summary?
#' @param show_group_counts Logical. Show N and event counts by group?
#' @param show_medians Logical. Show median survival by group?
#' @param show_ph_test Logical. Show proportional hazards test?
#' @param show_params Logical. Show original function parameters?
#' @param type Plot type: `"outcome"`, `"distribution"`, or `"forest"`.
#' @param reference_group Reference group for forest plot (e.g., `"G1"`).
#'
#' @examples
#' data(crc_virome)
#' res <- find_cutpoint(
#'   data = head(crc_virome, 50),
#'   predictor = "Alphapapillomavirus",
#'   outcome_time = "time_months",
#'   outcome_event = "status",
#'   num_cuts = 1,
#'   method = "systematic"
#' )
#'
#' @return An object of class `find_cutpoint` containing the
#'   optimal cut-points, statistic, and analysis parameters.
#'
#' @importFrom stats na.omit as.formula pchisq aggregate relevel
#' @importFrom survival Surv survfit survdiff coxph cox.zph
#' @importFrom cli cli_h1 cli_alert_info cli_alert_success cli_alert_warning
#' @importFrom cli cli_inform cli_abort cli_bullets cli_h2
#' @importFrom ggplot2 ggplot aes .data after_stat geom_density
#' @importFrom ggplot2 geom_histogram geom_vline labs theme_minimal
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom survminer ggsurvplot ggforest
#' @importFrom tools toTitleCase
#' @export
find_cutpoint <- function(data, predictor, outcome_time,
                          outcome_event,
                          num_cuts = 1, method = "systematic",
                          criterion = "logrank", covariates = NULL,
                          nmin = 20, seed = NULL, maxiter = 100,
                          quiet = FALSE, ...) {
  # --- 1. Validate User Inputs ---
  # Validates all user-facing arguments.
  .validate_find_cutpoint_inputs(
    data = data, predictor = predictor, outcome_time = outcome_time,
    outcome_event = outcome_event, num_cuts = num_cuts,
    method = method, criterion = criterion,
    covariates = covariates
  )

  # Store original name for the na_result helper
  original_predictor_name <- predictor

  # --- 2. Prepare Data ---
  # Subsets, removes NAs, and renames core columns
  # to "time", "event", "factor".
  userdata <- .prepare_cutpoint_data(
    data = data, predictor = predictor,
    outcome_time = outcome_time,
    outcome_event = outcome_event, covariates = covariates
  )

  # --- 3. Define the 'na_result' helper function ---
  # Standard object to return on failure.
  na_result <- function(userdata, num_cuts, method,
                        criterion, quiet) {
    output <- list(
      optimal_cuts = rep(NA, num_cuts),
      optimal_stat = NA,
      all_stats = NULL,
      userdata = userdata,
      parameters = list(
        method = method,
        analysis_type = "survival",
        predictor = original_predictor_name,
        num_cuts = num_cuts,
        criterion = criterion,
        covariates = covariates,
        nmin = nmin, # Store original nmin
        quiet = quiet
      )
    )
    class(output) <- "find_cutpoint"
    if (!quiet) {
      cli::cli_inform("Could not determine optimal cut-point.")
    }
    return(output)
  }

  # --- 4. Check for Empty Data (Post-NA removal) ---
  if (nrow(userdata) == 0) {
    if (!quiet) {
      cli::cli_inform("No complete cases found after removing NAs.")
    }
    return(na_result(
      userdata, num_cuts, method,
      criterion, quiet
    ))
  }

  # --- 5. Validate Data Conditions ---
  # Validates data: 0/1 status, non-negative time,
  # non-constant predictor, and sufficient sample size.
  validation_result <- .validate_data_conditions(
    userdata = userdata,
    nmin = nmin,
    num_cuts = num_cuts,
    outcome_event = outcome_event,
    quiet = quiet
  )

  if (!validation_result$valid) {
    return(na_result(
      userdata, num_cuts, method,
      criterion, quiet
    ))
  }

  # Get the calculated absolute nmin from the validation helper
  nmin_abs <- validation_result$nmin_abs

  # --- 6. Set Seed ---
  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  # --- 7. Route to appropriate search method ---
  output <- if (method == "systematic") {
    .systematic_search(
      userdata = userdata,
      num_cuts = num_cuts,
      criterion = criterion,
      covariates = covariates,
      nmin = nmin_abs, # Use calculated absolute value
      predictor_name = original_predictor_name,
      quiet = quiet,
      ...
    )
  } else { # genetic
    if (!quiet) {
      cli::cli_alert_info(paste(
        "Starting genetic search for",
        "{num_cuts} cut(s) using",
        "'{criterion}' criterion..."
      ))
    }

    ga_result <- .run_genetic_search(
      target = userdata$factor,
      numcut = num_cuts,
      time = userdata$time,
      censor = userdata$event,
      confound = if (!is.null(covariates)) {
        userdata[, covariates, drop = FALSE]
      } else {
        NULL
      },
      nmin = nmin_abs, # Use calculated absolute value
      criterion = criterion,
      numgen = maxiter,
      ...
    )

    # Check for the failure signal from .obj
    if (is.null(ga_result) ||
      !is.finite(ga_result$value) ||
      ga_result$value <= -.Machine$double.xmax) {
      if (!quiet) {
        cli::cli_inform("Genetic algorithm found no valid solution.")
      }
      return(na_result(
        userdata, num_cuts, method,
        criterion, quiet
      ))
    }

    optimal_cuts <- sort(ga_result$par[1:num_cuts])
    optimal_stat <- ga_result$value

    list(
      optimal_cuts = optimal_cuts,
      optimal_stat = optimal_stat,
      all_stats = NULL,
      userdata = userdata,
      parameters = list(
        method = "genetic",
        analysis_type = "survival",
        predictor = original_predictor_name,
        num_cuts = num_cuts,
        criterion = criterion,
        covariates = covariates,
        nmin = nmin, # Store original nmin
        quiet = quiet
      )
    )
  }

  class(output) <- "find_cutpoint"
  if (!quiet) print(output)
  invisible(output)
}

# --- Internal Helper: Systematic Search ---
#' @srrstats {RE1.0} Implements systematic grid search
#' for 1–2 cut-points.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#' @noRd
.systematic_search <- function(userdata, num_cuts, criterion,
                               covariates, nmin, predictor_name,
                               quiet, ...) {
  if (!quiet) {
    cli::cli_alert_info("Running systematic search...")
  }
  if (!quiet) cli::cli_alert_info("Testing for {num_cuts} cut-point(s)...")
  userdata <- userdata[order(userdata$factor), ]

  # Handle covariates
  cov_part <- if (!is.null(covariates)) {
    paste(" +", paste(covariates, collapse = " + "))
  } else {
    ""
  }

  # Pre-fit null model if criterion is "p_value"
  fit_null_model <- NULL
  if (criterion == "p_value") {
    # Fit model *without* factor to get baseline log-likelihood.
    null_formula_str <- if (!is.null(covariates)) {
      paste("Surv(time, event) ~", paste(covariates, collapse = " + "))
    } else {
      "Surv(time, event) ~ 1"
    }
    fit_null_model <- tryCatch(
      survival::coxph(as.formula(null_formula_str), data = userdata),
      error = function(e) NULL
    )
    if (is.null(fit_null_model)) {
      if (!quiet) {
        cli::cli_alert_warning(paste(
          "Could not fit null model for p-value. Aborting."
        ))
      }
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA,
        all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic",
          analysis_type = "survival",
          predictor = predictor_name,
          num_cuts = num_cuts,
          criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }
  }

  direction <- if (criterion == "p_value") "min" else "max"
  best_stat <- if (direction == "min") Inf else -Inf
  best_cut_val <- rep(NA, num_cuts)
  all_stats_df <- NULL

  if (num_cuts == 1) {
    search_grid <- unique(userdata$factor[nmin:(nrow(userdata) - nmin)])
    if (length(search_grid) == 0) {
      if (!quiet) {
        cli::cli_inform(paste(
          "Not enough data ({nrow(userdata)})",
          "for nmin ({nmin}) / {num_cuts}",
          "cut(s)."
        ))
      }
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA,
        all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic",
          analysis_type = "survival",
          predictor = predictor_name,
          num_cuts = num_cuts,
          criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }
    stats_per_cut <- vapply(search_grid, .get_stat,
      num_cuts = 1,
      data_in = userdata,
      criterion = criterion,
      cov_formula = cov_part, nmin = nmin,
      fit_null = fit_null_model, # Pass null model
      FUN.VALUE = numeric(1)
    )

    if (all(is.na(stats_per_cut))) {
      if (!quiet) {
        cli::cli_inform("No valid cut-points found (model failures).")
      }
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA,
        all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic",
          analysis_type = "survival",
          predictor = predictor_name,
          num_cuts = num_cuts,
          criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }
    best_idx <- if (direction == "min") {
      which.min(stats_per_cut)
    } else {
      which.max(stats_per_cut)
    }
    best_cut_val <- search_grid[best_idx]
    best_stat <- stats_per_cut[best_idx]
    all_stats_df <- data.frame(cut1 = search_grid, stat = stats_per_cut)
  } else {
    if (!quiet) {
      cli::cli_alert_info("Searching for 2 cuts is slow...")
    }

    # Set up sequential backend
    foreach::registerDoSEQ()

    possible_c1_indices <- nmin:(nrow(userdata) - (2 * nmin))
    grid1_values <- unique(userdata$factor[possible_c1_indices])

    if (length(grid1_values) == 0) {
      if (!quiet) {
        cli::cli_inform(paste(
          "Not enough data ({nrow(userdata)})",
          "for nmin ({nmin}) / {num_cuts}",
          "cut(s)."
        ))
      }
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA,
        all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic",
          analysis_type = "survival",
          predictor = predictor_name,
          num_cuts = num_cuts,
          criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }

    results_list <- foreach::foreach(
      c1 = grid1_values, .combine = "rbind",
      .export = c(".get_stat")
    ) %do% {
      best_local_stat <- if (direction == "min") Inf else -Inf
      best_local_c2 <- NA
      c1_max_index <- max(which(userdata$factor == c1))

      start_index_c2 <- c1_max_index + nmin
      end_index_c2 <- nrow(userdata) - nmin
      if (start_index_c2 > end_index_c2) {
        return(NULL)
      }
      possible_c2_indices <- start_index_c2:end_index_c2
      grid2_values <- unique(userdata$factor[possible_c2_indices])

      if (length(grid2_values) == 0) {
        return(NULL)
      }

      for (c2 in grid2_values) {
        stat <- .get_stat(c(c1, c2), 2, userdata, criterion,
          cov_part, nmin,
          fit_null = fit_null_model
        ) # Pass null model
        if (is.na(stat)) next

        is_better <- if (direction == "min") {
          (stat < best_local_stat)
        } else {
          (stat > best_local_stat)
        }
        if (is_better &&
          !is.infinite(stat)) {
          best_local_stat <- stat
          best_local_c2 <- c2
        }
      }
      if (is.na(best_local_c2)) {
        return(NULL)
      }
      data.frame(stat = best_local_stat, c1 = c1, c2 = c2)
    }

    if (is.null(results_list) || nrow(results_list) == 0) {
      if (!quiet) {
        cli::cli_inform(paste(
          "Not enough data ({nrow(userdata)})",
          "for nmin ({nmin}) / {num_cuts}",
          "cut(s)."
        ))
      }
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA,
        all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic",
          analysis_type = "survival",
          predictor = predictor_name,
          num_cuts = num_cuts,
          criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }

    best_idx <- if (direction == "min") {
      which.min(results_list$stat)
    } else {
      which.max(results_list$stat)
    }
    best_stat <- results_list$stat[best_idx]
    best_cut_val <- c(results_list$c1[best_idx], results_list$c2[best_idx])
  }

  output <- list(
    optimal_cuts = best_cut_val,
    optimal_stat = best_stat,
    all_stats = all_stats_df,
    userdata = userdata,
    parameters = list(
      method = "systematic",
      analysis_type = "survival",
      predictor = predictor_name,
      num_cuts = num_cuts,
      criterion = criterion,
      covariates = covariates,
      nmin = nmin,
      quiet = quiet
    )
  )
  if (!quiet) cli::cli_alert_success("Systematic search complete.")
  return(output)
}

# --- Internal Core Logic: Calculate Statistic ---
#' @srrstats {RE1.0} Uses `survdiff`/`coxph` for log-rank,
#' HR, p-value.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#' @noRd
.get_stat <- function(cuts, num_cuts, data_in,
                      criterion, cov_formula, nmin,
                      fit_null = NULL) {
  breaks <- sort(unique(c(-Inf, cuts, Inf)))
  num_intervals <- length(breaks) - 1
  data_in$group <- factor(cut(data_in$factor,
    breaks = breaks,
    labels = 1:num_intervals
  ))

  if (any(table(data_in$group) < nmin) ||
    nlevels(data_in$group) != (num_cuts + 1)) {
    return(NA)
  }

  formula_str <- paste("Surv(time, event) ~ group", cov_formula)

  if (criterion == "logrank") {
    # survdiff lacks covariate support; use Cox score test.
    if (cov_formula == "") {
      fit <- tryCatch(
        survival::survdiff(as.formula(formula_str), data = data_in),
        error = function(e) NULL
      )
      if (is.null(fit)) {
        return(NA)
      }
      return(fit$chisq)
    } else {
      # Use Cox score test (log-rank with covariates)
      fit <- tryCatch(
        survival::coxph(as.formula(formula_str), data = data_in),
        error = function(e) NULL
      )
      if (is.null(fit) ||
        is.null(fit$score)) {
        return(NA)
      }
      return(fit$score)
    }
  } else { # Cox-based criteria
    fit <- tryCatch(
      survival::coxph(as.formula(formula_str), data = data_in),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      return(NA)
    }

    if (criterion == "hazard_ratio") {
      if (is.null(fit$coefficients)) {
        return(NA)
      }
      # Get coef directly, don't use summary()
      coef_name <- paste0("group", num_cuts + 1)
      if (!coef_name %in% names(fit$coefficients)) {
        return(NA)
      }
      return(exp(fit$coefficients[coef_name]))
    } else if (criterion == "p_value") {
      # Use Likelihood Ratio Test
      if (is.null(fit_null) || is.null(fit$loglik)) {
        return(NA)
      }
      loglik0 <- fit_null$loglik[2] # Null model
      loglik1 <- fit$loglik[2] # Full model
      lrt_stat <- 2 * (loglik1 - loglik0)
      df_diff <- num_cuts # Df is number of groups - 1
      p_val <- stats::pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
      return(p_val)
    }
  }
}

# --- S3 Methods for Print, Summary, Plot ---
#' @rdname find_cutpoint
#' @srrstats {RE1.3, RE1.3a} PH diagnostics via `cox.zph`
#' in `summary()`.
#' @export
print.find_cutpoint <- function(x, ...) {
  if (is.null(x) ||
    any(is.na(x$optimal_cuts))) {
    cli::cli_inform("No optimal cut-point determined.")
    return(invisible(x))
  }

  method_name <- tools::toTitleCase(x$parameters$method)
  cli::cli_h1(paste(
    "Optimal Cut-point Analysis for Survival",
    "Data ({method_name})"
  ))

  stat_label <- switch(x$parameters$criterion,
    "logrank" = "Optimal Log-Rank Statistic",
    "hazard_ratio" = "Optimal Hazard Ratio",
    "p_value" = "Optimal P-value"
  )

  stat_val <- x$optimal_stat
  stat_val_fmt <- round(stat_val, 4)

  # For genetic p_value, the stat is the LRT statistic.
  if (x$parameters$criterion == "p_value" &&
    x$parameters$method == "genetic") {
    stat_label <- "Optimal LRT Statistic"
  }

  rounded_cuts <- round(x$optimal_cuts, 3)
  cli::cli_bullets(c(
    "*" = "Predictor: {.strong {x$parameters$predictor}}",
    "*" = "Criterion: {.strong {x$parameters$criterion}}",
    "*" = "{stat_label}: {.strong {stat_val_fmt}}",
    "v" = "Recommended Cut-point(s): {.strong {rounded_cuts}}"
  ))
  invisible(x)
}

#' @rdname find_cutpoint
#' @export
summary.find_cutpoint <- function(object,
                                  show_model = TRUE,
                                  show_group_counts = TRUE,
                                  show_medians = TRUE,
                                  show_ph_test = TRUE,
                                  show_params = TRUE, ...) {
  cli::cli_h1(paste(
    "Optimal Cut-point Analysis for Survival Data",
    "({tools::toTitleCase(object$parameters$method)})"
  ))

  if (is.null(object) ||
    any(is.na(object$optimal_cuts))) {
    cli::cli_inform("No valid optimal cut-point found.")
    return(invisible(object))
  }

  data <- object$userdata
  cuts <- object$optimal_cuts
  num_cuts <- object$parameters$num_cuts

  data$group <- factor(cut(data$factor,
    breaks = c(-Inf, cuts, Inf),
    labels = paste0("G", 1:(num_cuts + 1))
  ))

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
    fit_km <- survival::survfit(survival::Surv(time, event) ~ group,
      data = data
    )
    print(fit_km)
  }

  formula_str <- "survival::Surv(time, event) ~ group"
  if (!is.null(object$parameters$covariates)) {
    formula_str <- paste(
      formula_str, "+",
      paste(object$parameters$covariates,
        collapse = " + "
      )
    )
  }
  fit_cox <- tryCatch(
    survival::coxph(as.formula(formula_str), data = data),
    error = function(e) NULL
  )

  if (show_model) {
    cli::cli_h2("Final Cox Model Summary")
    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model: convergence failed.")
    } else {
      print(summary(fit_cox))
    }
  }

  if (show_ph_test && !is.null(fit_cox)) {
    cli::cli_h2("Proportional Hazards Assumption Test")
    print(survival::cox.zph(fit_cox))
  }

  if (show_params) {
    cli::cli_h1("Analysis Parameters")
    params <- object$parameters
    param_bullets <- c(
      "*" = "Search Method: {tools::toTitleCase(params$method)}",
      "*" = "Predictor: {params$predictor}",
      "*" = "Number of cuts: {params$num_cuts}",
      "*" = "Minimum group size (nmin): {params$nmin}"
    )
    if (!is.null(params$covariates)) {
      param_bullets <- c(
        param_bullets,
        "*" = "Covariates: {paste(params$covariates, collapse = ', ')}"
      )
    }
    cli::cli_bullets(param_bullets)
  }

  invisible(object)
}

#' @rdname find_cutpoint
#' @export
plot.find_cutpoint <- function(x, type = "outcome",
                               reference_group = NULL, ...) {
  type <- match.arg(type, choices = c("outcome", "distribution", "forest"))

  if (is.null(x) || any(is.na(x$optimal_cuts))) {
    cli::cli_inform("Cannot generate plot: no valid cut-point found.")
    return(invisible(NULL))
  }

  cuts <- x$optimal_cuts
  data <- x$userdata
  num_cuts <- length(cuts)
  data$group <- factor(cut(data$factor,
    breaks = c(-Inf, cuts, Inf),
    labels = paste0("G", 1:(num_cuts + 1))
  ))

  if (type == "distribution") {
    p <- ggplot2::ggplot(x$userdata, ggplot2::aes(x = .data$factor)) +
      ggplot2::geom_histogram(
        aes(y = ggplot2::after_stat(density)),
        bins = 30,
        fill = "#56B4E9", color = "black", alpha = 0.7
      ) +
      ggplot2::geom_density(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_vline(
        xintercept = cuts, color = "#D55E00",
        linetype = "dashed", linewidth = 1.2
      ) +
      ggplot2::labs(
        title = "Distribution of Predictor with Optimal Cut-points",
        x = x$parameters$predictor, y = "Density"
      ) +
      ggplot2::theme_minimal()
    return(p)
  } else if (type == "outcome") {
    fit <- survival::survfit(Surv(time, event) ~ group, data = data)
    p <- survminer::ggsurvplot(fit,
      data = data,
      pval = TRUE,
      risk.table = TRUE,
      legend.title = "Groups",
      palette = "jco",
      ggtheme = ggplot2::theme_minimal()
    )
    p$plot <- p$plot + ggplot2::labs(
      title = paste("Survival Curves by", x$parameters$predictor, "Group")
    )
    return(p)
  } else if (type == "forest") {
    if (!requireNamespace("broom", quietly = TRUE)) {
      stop("Package 'broom' is required for the forest plot.",
        call. = FALSE
      )
    }

    cli::cli_alert_info("Generating Forest Plot of Hazard Ratios...")

    group_levels <- levels(data$group)
    if (is.null(reference_group) || !(reference_group %in% group_levels)) {
      reference_group <- group_levels[1]
      cli::cli_alert_warning(
        paste(
          "Invalid/missing reference group.",
          "Defaulting to: {reference_group}"
        )
      )
    }
    data$group <- stats::relevel(data$group, ref = reference_group)

    formula_str <- "survival::Surv(time, event) ~ group"
    if (!is.null(x$parameters$covariates)) {
      formula_str <- paste(
        formula_str, "+",
        paste(x$parameters$covariates,
          collapse = " + "
        )
      )
    }
    fit_cox <- tryCatch(
      survival::coxph(as.formula(formula_str), data = data),
      error = function(e) NULL
    )

    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model for forest plot.")
      return(invisible(NULL))
    }

    p <- survminer::ggforest(fit_cox, data = data) +
      ggplot2::labs(
        title = "Hazard Ratios for Predictor Groups",
        subtitle = paste("Reference group:", reference_group),
        x = "Hazard Ratio (95% CI)"
      )
    return(p)
  }
}
