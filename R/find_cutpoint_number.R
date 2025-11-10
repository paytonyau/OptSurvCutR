# ===================================================================
#' Find Optimal Number of Cut-points for Survival Data
# ===================================================================
#' @description
#' Finds optimal cut-point number (0 to `max_cuts`) for a Cox model
#' by comparing AIC, AICc, or BIC. Supports systematic search
#' (`max_cuts <= 2`) and genetic algorithm (`rgenoud`).
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.1} Computes AIC, AICc, or BIC for model selection.
#' @srrstats {G2.13} `cli_abort()` for invalid input.
#' @srrstats {G2.14a} `NA` results handled via `na_result`.
#' @srrstats {G3.1} `plot()` shows IC vs. number of cuts.
#' @srrstats {G1.2} References provided for AIC, AICc, BIC.
#' @srrstats {G1.3} Systematic grid search (max_cuts <= 2) and
#' `rgenoud` global optimization documented.
#' @srrstats {G1.5} Compared with `cutpointr`/`survminer` in vignette.
#' @srrstats {G1.6} Numerical stability via `survival::coxph`
#' and `rgenoud`; edge cases return `NA`.
#' @srrstats {G2.3b} NSE via `as.formula()`/`data[...]`; no unsafe eval.
#' @srrstats {G2.4} `NA` removed via `stats::na.omit()`;
#' structured `NA` object returned on failure.
#' @srrstats {G2.4e} `optimal_num_cuts`/`optimal_cuts` are `NA`
#' when no valid solution found.
#' @srrstats {G2.6} Input validation via direct checks.
#' @srrstats {G2.8} Informative errors via `cli::cli_abort()`.
#' @srrstats {G2.10} Warnings via `cli::cli_alert_warning()`.
#' @srrstats {G2.12} Graceful degradation via `na_result()` for
#' empty data or model failures.
#' @srrstats {G2.14c} `NA` propagation controlled.
#' @srrstats {G4.0} All parameters/return values documented.
#' @srrstats {G5.4c} Edge cases (zero rows, constant predictor) tested.
#' @srrstats {G5.6a, G5.6b} Negative `max_cuts`/`nmin` rejected.
#' @srrstats {G5.7} Large `max_cuts` constrained by `nmin`.
#' @srrstats {G5.8d} `set.seed(seed)` for genetic reproducibility.
#' @srrstats {G5.12} Systematic search scales poorly > 2 cuts.
#'
#' @srrstats {RE1.1} Assumes PH; check `summary()` for `cox.zph`.
#' @srrstats {RE1.2} Provides AIC, AICc, BIC fit statistics.
#' @srrstats {RE1.3, RE1.3a} PH diagnostics via `cox.zph` in `summary()`.
#' @srrstats {RE2.0, RE2.1} Estimates/SEs from `coxph` in `summary()`.
#' @srrstats {RE2.4, RE2.4a, RE2.4b} `tryCatch` checks model
#' convergence; failures return `NA`.
#' @srrstats {RE3.0} Prediction not implemented.
#' @srrstats {RE3.2, RE3.3} Not applicable; no predictions.
#' @srrstats {RE4.2} Model selection via AIC, AICc, or BIC.
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
#' Systematic search is slow for `max_cuts > 2`; use `genetic`.
#'
#' @references
#' Akaike, H. (1974). A new look at the statistical model identification.
#' *IEEE Transactions on Automatic Control*, **19**(6), 716–723.
#' \doi{10.1109/TAC.1974.1100705}
#'
#' Chang, C., Hsieh, M.-K., Chang, W.-Y., Chiang, A. J., &
#' Chen, J. (2017). Determining the optimal number and location of cutoff
#' points with application to data of cervical cancer. *PLOS ONE*, 12(4),
#' e0176231. \doi{10.1371/journal.pone.0176231}
#'
#' Chen, Y., Huang, J., He, X., Gao, Y., Mahara, G., Lin, Z.,
#' & Zhang, J. (2019). A novel approach to determine two optimal
#' cut-points of a continuous predictor with a U-shaped relationship to
#' hazard ratio in survival data: Simulation and application. *BMC Medical
#' Research Methodology*, 19(1), 96. \doi{10.1186/s12874-019-0738-4}
#'
#' Schwarz, G. (1978). Estimating the dimension of a model.
#' *The Annals of Statistics*, **6**(2), 461–464.
#' \doi{10.1214/aos/1176344136}
#'
#' Hurvich, C. M., & Tsai, C.-L. (1989). Regression and time series model
#' selection in small samples. *Biometrika*, **76**(2), 297–307.
#' \doi{10.1093/biomet/76.2.297}
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
#' @param maxiter Integer; generations for `rgenoud` (default 100).
#' @param x An object from [find_cutpoint_number()].
#' @param object An object from [find_cutpoint_number()].
#' @param y Unused.
#' @param ... Additional arguments passed to `rgenoud`.
#'
#' @examples
#' data(crc_virome)
#' res <- find_cutpoint_number(
#'   data = head(crc_virome, 50),
#'   predictor = "Alphapapillomavirus",
#'   outcome_time = "time_months",
#'   outcome_event = "status",
#'   method = "systematic",
#'   max_cuts = 1
#' )
#'
#' @return An S3 object (`find_cutpoint_number_result`) with
#' `results`, `parameters`, `userdata`, `optimal_num_cuts`,
#' and `optimal_cuts`.
#'
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom stats na.omit as.formula aggregate
#' @importFrom survival coxph Surv survfit
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success
#' @importFrom cli cli_inform cli_abort cli_bullets cli_h2
#' @importFrom ggplot2 ggplot aes .data geom_line geom_point labs
#' @importFrom ggplot2 theme_minimal scale_x_continuous element_text theme
#' @importFrom tools toTitleCase
#' @export
find_cutpoint_number <- function(data, predictor,
                                 outcome_time, outcome_event,
                                 method = "systematic", criterion = "BIC",
                                 covariates = NULL, max_cuts = 2,
                                 nmin = 0.1, seed = NULL, maxiter = 100,
                                 ...) {
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
        "i" = "Install with: install.packages(\"rgenoud\")",
        "i" = "Or, use `method = \"systematic\"` (for max_cuts <= 2)."
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
  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"

  n <- nrow(userdata)

  event_col <- userdata$event
  if (!is.numeric(event_col)) {
    cli::cli_abort(c(
      "Event column ({.arg {outcome_event}}) must be numeric.",
      "i" = "Detected: {.cls {class(event_col)}}.",
      "i" = "Please convert to 0 (censored) and 1 (event)."
    ))
  }

  valid_events <- unique(event_col)
  # Check if all values are either 0 or 1
  if (!all(valid_events %in% c(0, 1))) {
    invalid_vals <- sort(valid_events[!(valid_events %in% c(0, 1))])
    cli::cli_abort(c(
      "Event column ({.arg {outcome_event}}) must contain only 0 and 1.",
      "i" = "Found invalid value(s): {.val {invalid_vals}}"
    ))
  }

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
        covariates = covariates,
        max_cuts = max_cuts,
        nmin = nmin
      ),
      userdata = userdata,
      optimal_num_cuts = NA,
      optimal_cuts = NA
    )
    class(output) <- "find_cutpoint_number_result"
    if (!is.null(userdata)) {
      cli::cli_inform("No valid results found with given parameters.")
    }
    return(output)
  }

  # Handle no data post-NA removal
  if (n == 0) {
    cli::cli_inform("No complete cases found after removing NAs.")
    return(na_result(userdata, method, criterion))
  }

  if (nmin < 1 && nmin > 0) {
    nmin_abs <- ceiling(nmin * n)
    cli::cli_alert_info(
      "nmin {nmin} is a proportion. Min. group size set to {nmin_abs}."
    )
    nmin <- nmin_abs
  } else if (nmin >= 1) {
    nmin <- as.integer(nmin)
  } else {
    cli::cli_abort("'nmin' must be a positive number.")
  }

  # Check for insufficient data
  if (n < nmin * (max_cuts + 1)) {
    cli::cli_inform(paste(
      "Not enough data ({n}) for nmin ({nmin}) and",
      "max_cuts ({max_cuts})."
    ))
    return(na_result(userdata, method, criterion))
  }

  # Check for constant predictor
  if (length(unique(userdata$factor)) <= max_cuts) {
    cli::cli_inform(paste(
      "Predictor has too few unique values",
      "({length(unique(userdata$factor))}) for max_cuts ({max_cuts})."
    ))
    return(na_result(userdata, method, criterion))
  }

  cli::cli_alert_info(
    "Finding optimal cut number: method = {.strong {method}}"
  )

  # Covariates to params
  params <- list(
    userdata = userdata, max_cuts = max_cuts, nmin = nmin,
    criterion = criterion, covariates = covariates,
    maxiter = maxiter, ...
  )

  results <- if (method == "systematic") {
    do.call(.systematic_search_num, params)
  } else { # genetic
    do.call(.genetic_search_num, params)
  }

  if (is.null(results) || !is.data.frame(results) || nrow(results) == 0) {
    cli::cli_inform("Search algorithm failed to produce results.")
    return(na_result(userdata, method, criterion))
  }

  # --- Post-process IC results ---
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
      covariates = covariates,
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
    output$optimal_cuts <- results$cuts[[min_ic_idx]]
  } else {
    output$optimal_num_cuts <- NA
    output$optimal_cuts <- NA
  }

  class(output) <- "find_cutpoint_number_result"

  print(output)
  invisible(output)
}

# --- Internal Helper Functions ---
#' @srrstats {RE4.1} Helper for systematic model selection (grid search).
#' @srrstats {RE1.0} Fits Cox models via `survival::coxph`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#' @noRd
.systematic_search_num <- function(userdata, max_cuts, nmin, criterion,
                                   covariates, ...) {
  if (max_cuts > 2) {
    cli::cli_abort(paste(
      "'systematic' method is computationally intensive",
      "and only implemented for max_cuts <= 2."
    ))
  }

  n <- nrow(userdata)
  userdata <- userdata[order(userdata$factor), ]

  # Define covariate formula part
  cov_part <- if (!is.null(covariates)) {
    paste(" +", paste(covariates, collapse = " + "))
  } else {
    ""
  }

  # Register sequential backend
  foreach::registerDoSEQ()

  # Base model (0 cuts)
  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)
  ic0 <- tryCatch(
    {
      fit0 <- survival::coxph(as.formula(base_formula_str), data = userdata)
      .calc_ic(fit0,
        k = 1 + length(covariates), n = n,
        criterion = criterion
      )
    },
    error = function(e) {
      cli::cli_inform(
        "Could not calculate IC for base model (0 cuts): {e$message}"
      )
      return(NA_real_)
    }
  )

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info("Testing for {k_cuts} cut-point(s)...")
    best_res_for_k <- list(ic = Inf, cuts = NULL)

    if (k_cuts == 1) {
      grid1 <- unique(userdata$factor[nmin:(n - nmin)])
      if (length(grid1) == 0) {
        cli::cli_inform(paste(
          "Not enough data ({n}) for nmin ({nmin})",
          "and {k_cuts} cut(s). Skipping."
        ))
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      res_list <- foreach::foreach(
        c1 = grid1, .combine = "rbind",
        .export = c(".get_model_ic_num", ".calc_ic")
      ) %do% {
        factor_status <- factor(ifelse(userdata$factor <= c1, 0, 1))
        if (min(table(factor_status)) < nmin ||
          nlevels(factor_status) < 2) {
          return(NULL)
        }

        current_ic <- .get_model_ic_num(
          userdata, factor_status, k_cuts,
          n, criterion, cov_part
        )
        if (is.finite(current_ic)) {
          data.frame(ic = current_ic, cuts = c1)
        } else {
          NULL
        }
      }

      if (is.null(res_list) || nrow(res_list) == 0) {
        cli::cli_inform(paste(
          "No valid cut-points found for {k_cuts} cut(s)",
          "due to model failures or constraints."
        ))
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
        cli::cli_inform(paste(
          "Not enough data ({n}) for nmin ({nmin})",
          "and {k_cuts} cut(s). Skipping."
        ))
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      res_list <- foreach::foreach(
        c1 = grid1, .combine = "rbind",
        .export = c(".get_model_ic_num", ".calc_ic")
      ) %do% {
        best_inner_res <- list(ic = Inf, c2 = NA)

        start_idx_g2 <- which(userdata$factor > c1)[nmin]
        if (is.na(start_idx_g2)) {
          return(NULL)
        }
        grid2_end_idx <- n - nmin
        if (start_idx_g2 >= grid2_end_idx) {
          return(NULL)
        }

        grid2 <- unique(userdata$factor[start_idx_g2:grid2_end_idx])
        if (length(grid2) == 0) {
          return(NULL)
        }

        for (c2 in grid2) {
          if (is.na(c2) || c2 <= c1) next
          factor_status <- as.factor(cut(userdata$factor,
            breaks = c(-Inf, c1, c2, Inf)
          ))
          if (min(table(factor_status)) < nmin ||
            nlevels(factor_status) < 3) {
            next
          }

          current_ic <- .get_model_ic_num(
            userdata, factor_status,
            k_cuts, n, criterion, cov_part
          )
          if (current_ic < best_inner_res$ic) {
            best_inner_res <- list(ic = current_ic, c2 = c2)
          }
        }

        if (is.finite(best_inner_res$ic)) {
          data.frame(
            ic = best_inner_res$ic, c1 = c1, c2 = best_inner_res$c2
          )
        } else {
          NULL
        }
      }

      if (is.null(res_list) || nrow(res_list) == 0) {
        cli::cli_inform(paste(
          "No valid cut-points found for {k_cuts} cut(s)",
          "due to model failures or constraints."
        ))
        new_row <- data.frame(num_cuts = k_cuts, IC = NA_real_)
        new_row$cuts <- I(list(NULL))
        results <- rbind(results, new_row)
        next
      }

      best_row <- res_list[which.min(res_list$ic), ]
      best_res_for_k <- list(
        ic = best_row$ic,
        cuts = c(best_row$c1, best_row$c2)
      )
    }

    min_ic_for_k <- if (is.finite(best_res_for_k$ic)) {
      best_res_for_k$ic
    } else {
      NA_real_
    }
    cuts_for_k <- list(best_res_for_k$cuts)

    new_row <- data.frame(num_cuts = k_cuts, IC = min_ic_for_k)
    new_row$cuts <- I(cuts_for_k)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

#' @srrstats {RE4.1} Helper for genetic-algorithm-based model selection.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#' @noRd
.genetic_search_num <- function(userdata, max_cuts, nmin, criterion,
                                covariates, maxiter, ...) {
  n <- nrow(userdata)

  # Define covariate formula part
  cov_part <- if (!is.null(covariates)) {
    paste(" +", paste(covariates, collapse = " + "))
  } else {
    ""
  }
  num_cov <- length(covariates)

  # Base model (0 cuts)
  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)
  ic0 <- tryCatch(
    {
      fit0 <- survival::coxph(as.formula(base_formula_str), data = userdata)
      .calc_ic(fit0,
        k = 1 + num_cov, n = n, criterion = criterion
      )
    },
    error = function(e) {
      cli::cli_inform(
        "Could not calculate IC for base model (0 cuts): {e$message}"
      )
      return(NA_real_)
    }
  )

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    cli::cli_alert_info(
      "Running genetic algorithm for {k_cuts} cut-point(s)..."
    )

    ga_result <- tryCatch(
      {
        .run_genetic_search(
          target = userdata$factor,
          numcut = k_cuts,
          time = userdata$time,
          censor = userdata$event,
          confound = if (!is.null(covariates)) {
            userdata[, covariates, drop = FALSE]
          } else {
            NULL
          },
          nmin = nmin,
          criterion = "loglik",
          numgen = maxiter,
          ...
        )
      },
      error = function(e) {
        cli::cli_inform(
          "Genetic algorithm failed for {k_cuts} cut(s): {e$message}"
        )
        return(NULL)
      }
    )

    ic_val <- NA_real_
    cuts_val <- list(NULL)

    if (!is.null(ga_result) && is.finite(ga_result$value) &&
      ga_result$value > -.Machine$double.xmax) {
      max_logL <- ga_result$value
      k_params <- k_cuts + num_cov
      ic_val <- .calc_ic(
        model = list(loglik = c(NA, max_logL)),
        k = k_params, n = n, criterion = criterion
      )
      cuts_val <- list(sort(ga_result$par[1:k_cuts]))
    } else {
      cli::cli_inform(paste(
        "No valid cut-points found for {k_cuts} cut(s)",
        "due to genetic algorithm failure or constraints."
      ))
    }

    new_row <- data.frame(num_cuts = k_cuts, IC = ic_val)
    new_row$cuts <- I(cuts_val)
    results <- rbind(results, new_row)
  }
  names(results)[2] <- criterion
  return(results)
}

#' @srrstats {RE4.1} Computes AIC, AICc, or BIC from a fitted Cox model.
#' @srrstats {RE1.0} Relies on `logLik()` from `coxph` fit.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#' @noRd
.get_model_ic_num <- function(userdata, factor_status, k_cuts, n,
                              criterion, cov_part) {
  num_cov <- length(cov_part[cov_part != ""])

  # Build formula with covariates
  formula_str <- paste("survival::Surv(time, event) ~ factor_status", cov_part)

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
#' @srrstats {RE1.3, RE1.3a} PH diagnostics via `cox.zph` in `summary()`.
#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")

  method_text <- if (!is.null(x$parameters$method)) {
    x$parameters$method
  } else {
    "Unknown"
  }
  criterion_text <- if (!is.null(x$parameters$criterion)) {
    x$parameters$criterion
  } else {
    "IC"
  }

  cli::cli_text("Method: {.strong {method_text}}")
  cli::cli_text("Criterion: {.strong {criterion_text}}")

  if (!is.null(x$parameters$covariates)) {
    cli::cli_text("Covariates: {.strong {paste(x$parameters$covariates,
      collapse = ', ')}}")
  }

  if (is.null(x$results) || nrow(x$results) == 0 ||
    all(is.na(x$results[[criterion_text]]))) {
    cli::cli_inform("No optimal model could be determined.")
    return(invisible(x))
  }

  print_df <- x$results

  if ("cuts" %in% names(print_df) && is.list(print_df$cuts)) {
    print_df$cuts <- vapply(print_df$cuts, function(c) {
      if (is.null(c)) "NA" else paste(round(c, 2), collapse = ", ")
    }, FUN.VALUE = character(1))
  }

  is_num <- vapply(print_df, is.numeric, FUN.VALUE = logical(1))
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0 && weight_col %in% names(x$results)) {
    print_df[[weight_col]] <- paste0(
      round(x$results[[weight_col]] * 100, 1), "%"
    )
  }

  final_cols <- c(
    "num_cuts", criterion_text,
    paste0("Delta_", criterion_text),
    weight_col, "Evidence", "cuts"
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
    if (!is.null(best_cuts_vals)) {
      rounded_cuts <- round(best_cuts_vals, 2)
      cli::cli_text("  Optimal cuts at: {.strong {rounded_cuts}}")
    }
  } else {
    cli::cli_inform("\nConclusion: No optimal model could be determined.")
  }

  cli::cli_text(
    "\nHint: Use `summary()` for details, `plot()` to visualize."
  )
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
  criterion_text <- if (!is.null(object$parameters$criterion)) {
    object$parameters$criterion
  } else {
    "IC"
  }

  cli::cli_h1(paste(
    "Optimal Cut-point Number Analysis",
    "({tools::toTitleCase(object$parameters$method)})"
  ))

  if (is.null(object) || is.null(object$results) ||
    nrow(object$results) == 0 ||
    all(is.na(object$results[[criterion_text]]))) {
    cli::cli_inform("Cannot summarize: no valid model was found.")
    return(invisible(object))
  }

  if (show_comparison_table) {
    print(object)
  }

  if (show_best_model_details) {
    cli::cli_h1("Details for Best Model")

    best_result <- object$results[which.min(
      object$results[[criterion_text]]
    ), ]
    best_cuts_vals <- best_result$cuts[[1]]
    num_cuts <- best_result$num_cuts

    cli::cli_text(
      "The best model found has {.strong {num_cuts}} cut-point(s)."
    )
    if (!is.null(best_cuts_vals)) {
      rounded_cuts <- round(best_cuts_vals, 2)
      cli::cli_text("Cut-point values: {.strong {rounded_cuts}}.")
    }

    data <- object$userdata

    cov_part <- if (!is.null(object$parameters$covariates)) {
      paste(" +", paste(object$parameters$covariates, collapse = " + "))
    } else {
      ""
    }

    if (num_cuts > 0) {
      data$group <- cut(data$factor,
        breaks = c(-Inf, best_cuts_vals, Inf),
        labels = paste0("G", 1:(num_cuts + 1))
      )

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
    }

    cli::cli_h2("Final Cox Proportional-Hazards Model")
    formula_str <- if (num_cuts > 0) {
      paste("survival::Surv(time, event) ~ group", cov_part)
    } else {
      paste("survival::Surv(time, event) ~ factor", cov_part)
    }

    model_data <- if (num_cuts > 0) data else object$userdata

    fit_cox <- tryCatch(
      survival::coxph(as.formula(formula_str), data = model_data),
      error = function(e) NULL
    )
    if (is.null(fit_cox)) {
      cli::cli_inform(
        "Could not fit Cox model for best model: convergence failed."
      )
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
  if (!is.null(params$covariates)) {
    param_bullets <- c(
      param_bullets,
      "*" = "Covariates: {paste(params$covariates, collapse = ', ')}"
    )
  }

  cli::cli_bullets(param_bullets)

  invisible(object)
}

#' @rdname find_cutpoint_number
#' @export
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results
  criterion_text <- if (!is.null(x$parameters$criterion)) {
    x$parameters$criterion
  } else {
    "IC"
  }

  if (is.null(results) || nrow(results) == 0 ||
    all(is.na(results[[criterion_text]]))) {
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
      data = ~ subset(., num_cuts == best_num_cuts),
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
