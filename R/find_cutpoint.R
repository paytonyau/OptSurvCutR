#' Find Optimal Cut-points for Survival Data
#'
#' @description
#' Determines optimal cut-point(s) for a continuous predictor in a time-to-event
#' (survival) analysis context. It can use a systematic search for 1-2 cut-points
#' or a more flexible genetic algorithm for any number of cut-points.
#'
#' @param data A data frame containing the analysis variables.
#' @param predictor The name of the continuous predictor variable to find cut-points for.
#' @param outcome_time The name of the time-to-event variable for survival analysis.
#' @param outcome_event The name of the event status variable for survival analysis.
#' @param num_cuts The number of cut-points to find. Default is 1.
#' @param method The algorithm to use: "systematic" or "genetic".
#' @param criterion The statistical criterion to optimize. Options are:
#'   "logrank" (maximizes the log-rank chi-squared statistic),
#'   "hazard_ratio" (maximizes the Hazard Ratio from a Cox model), or
#'   "p_value" (minimizes the p-value from a Cox model).
#' @param covariates A character vector of covariate names to include in the model.
#' @param nmin The minimum number of observations required in each group created by the cut-points.
#'   Can be specified as an integer (e.g., 20) or a proportion (e.g., 0.1).
#' @param seed An optional integer for setting the random seed to ensure
#'   reproducible results when using the "genetic" method.
#' @param maxiter The number of generations for the genetic algorithm. Default is 100.
#' @param quiet Logical. If TRUE, suppresses the printing of the final result object.
#' @param ... Additional arguments passed to the genetic algorithm (e.g., `popSize`).
#' @param x An object from \code{\link{find_cutpoint}}.
#' @param object An object from \code{\link{find_cutpoint}}.
#' @param show_model Logical. If TRUE, shows the full summary of the final Cox model.
#' @param show_group_counts Logical. If TRUE, shows the number of subjects and events in each group.
#' @param show_medians Logical. If TRUE, shows the median survival for each group.
#' @param show_ph_test Logical. If TRUE, shows the proportional hazards assumption test.
#' @param show_params Logical. If TRUE, shows the parameters of the original function call.
#' @param type The type of plot to generate: "outcome" (a survival plot), "distribution", or "forest".
#' @param reference_group The reference group for the forest plot (e.g., "G1").
#'
#' @return An object of class \code{find_cutpoint} containing the optimal cut-points,
#'   the corresponding statistic, and other parameters used in the analysis.
#' @importFrom stats na.omit as.formula pchisq sd rnorm runif anova aggregate relevel
#' @importFrom survival Surv survfit survdiff coxph cox.zph
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_alert_danger cli_alert_warning cli_inform cli_abort
#' @importFrom ggplot2 ggplot aes .data geom_line geom_vline labs theme_minimal
#'   geom_histogram geom_density
#' @importFrom foreach %do% registerDoSEQ
#' @importFrom survminer ggsurvplot ggforest
#' @importFrom tools toTitleCase
#' @export
find_cutpoint <- function(data, predictor, outcome_time, outcome_event, num_cuts = 1, method = "systematic",
                          criterion = "logrank", covariates = NULL, nmin = 20, seed = NULL, maxiter = 100, quiet = FALSE, ...) {

  # --- 1. Input Validation and Data Prep ---
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("logrank", "hazard_ratio", "p_value"))

  if (!is.numeric(num_cuts) || num_cuts < 0 || num_cuts != round(num_cuts)) {
    stop("num_cuts must be a non-negative integer", call. = FALSE)
  }
  if (criterion == "hazard_ratio" && num_cuts > 1) {
    stop("'hazard_ratio' criterion is only supported for num_cuts = 1.", call. = FALSE)
  }

  # Graceful failure if rgenoud is missing
  if (method == "genetic" && !requireNamespace("rgenoud", quietly = TRUE)) {
    cli::cli_abort(
      c("The 'genetic' method requires the 'rgenoud' package.",
        "i" = "Please install it by running: install.packages(\"rgenoud\")",
        "i" = "Alternatively, use `method = \"systematic\"` (for num_cuts <= 2).")
    )
  }

  # Check for systematic search constraints
  if (method == "systematic" && !num_cuts %in% c(1, 2)) {
    stop("Systematic search currently only supports num_cuts = 1 or 2.", call. = FALSE)
  }

  # Validate required columns
  if (is.null(predictor)) stop("A 'predictor' variable must be specified.", call. = FALSE)
  if (is.null(outcome_time) || is.null(outcome_event)) stop("Both 'outcome_time' and 'outcome_event' must be specified.", call. = FALSE)

  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  missing_vars <- setdiff(required_vars, names(data))
  if (length(missing_vars) > 0) {
    stop(paste0("The following required column(s) were not found in the data: '",
                paste(missing_vars, collapse = "', '"), "'"), call. = FALSE)
  }

  # Data preparation
  userdata <- data[, unique(required_vars), drop = FALSE]
  userdata <- stats::na.omit(userdata)

  # --- Create a standard NA result object ---
  na_result <- function(userdata, num_cuts, method, criterion, quiet) {
    output <- list(
      optimal_cuts = rep(NA, num_cuts),
      optimal_stat = NA,
      all_stats = NULL,
      userdata = userdata,
      parameters = list(
        method = method,
        analysis_type = "survival",
        predictor = predictor,
        num_cuts = num_cuts,
        criterion = criterion,
        covariates = covariates,
        nmin = nmin,
        quiet = quiet
      )
    )
    class(output) <- "find_cutpoint"
    if (!quiet) cli::cli_inform("No optimal cut-point could be determined with the given parameters.")
    return(output)
  }

  # Handle cases with no data after NA removal
  if (nrow(userdata) == 0) {
    if (!quiet) cli::cli_inform("No complete cases found in the data after removing NAs.")
    return(na_result(userdata, num_cuts, method, criterion, quiet))
  }

  original_predictor_name <- predictor
  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"

  if (any(userdata$time < 0, na.rm = TRUE)) {
    stop("Time variable must be non-negative.", call. = FALSE)
  }

  # Handle nmin as proportion or integer
  if (nmin > 0 && nmin < 1) {
    nmin_abs <- floor(nmin * nrow(userdata))
    if (!quiet) cli::cli_alert_info("Interpreting nmin = {nmin} as a proportion. Minimum group size set to {nmin_abs}.")
    nmin <- nmin_abs
  } else if (nmin < 0) {
    stop("'nmin' must be a non-negative number.", call. = FALSE)
  }

  # Check for insufficient data
  if (nrow(userdata) < nmin * (num_cuts + 1)) {
    if (!quiet) cli::cli_inform("Not enough data ({nrow(userdata)}) for nmin ({nmin}) and {num_cuts} cut(s). Returning NA cuts.")
    return(na_result(userdata, num_cuts, method, criterion, quiet))
  }

  # Check for constant predictor
  if (length(unique(userdata$factor)) <= num_cuts) {
    if (!quiet) cli::cli_inform("Predictor has too few unique values ({length(unique(userdata$factor))}) for {num_cuts} cut(s). Returning NA cuts.")
    return(na_result(userdata, num_cuts, method, criterion, quiet))
  }

  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  # --- 2. Route to appropriate search method ---
  output <- if (method == "systematic") {
    .systematic_search(
      userdata = userdata,
      num_cuts = num_cuts,
      criterion = criterion,
      covariates = covariates,
      nmin = nmin,
      predictor_name = original_predictor_name,
      quiet = quiet,
      ...
    )
  } else { # genetic
    if (!quiet) cli::cli_alert_info("Starting genetic algorithm for {num_cuts} cut-point(s) using '{criterion}' criterion...")

    ga_result <- .run_genetic_search(
      target = userdata$factor,
      numcut = num_cuts,
      time = userdata$time,
      censor = userdata$event,
      confound = if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL,
      nmin = nmin,
      criterion = criterion,
      numgen = maxiter,
      ...
    )

    if (is.null(ga_result) || !is.finite(ga_result$value)) {
      if (!quiet) cli::cli_inform("Genetic algorithm could not find a valid solution with the given constraints.")
      return(na_result(userdata, num_cuts, method, criterion, quiet))
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
        nmin = nmin,
        quiet = quiet
      )
    )
  }

  class(output) <- "find_cutpoint"
  if (!quiet) print(output)
  invisible(output)
}

# --- Internal Helper: Systematic Search ---
.systematic_search <- function(userdata, num_cuts, criterion, covariates, nmin, predictor_name, quiet, ...) {
  if (!quiet) cli::cli_alert_info("Finding optimal number of cuts: method = systematic")
  if (!quiet) cli::cli_alert_info("Testing for {num_cuts} cut-point(s)...")
  userdata <- userdata[order(userdata$factor), ]

  # Handle covariates
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

  # Pre-fit null model if criterion is "p_value"
  fit_null_model <- NULL
  if (criterion == "p_value") {
    # We fit a model *without* the factor to get the baseline log-likelihood.
    # If there are no covariates, the null model is intercept-only.
    null_formula_str <- if (!is.null(covariates)) paste("Surv(time, event) ~", paste(covariates, collapse = " + ")) else "Surv(time, event) ~ 1"
    fit_null_model <- tryCatch(survival::coxph(as.formula(null_formula_str), data = userdata), error = function(e) NULL)
    if (is.null(fit_null_model)) {
      if (!quiet) cli::cli_alert_warning("Could not fit null model for p-value calculation. Aborting systematic search.")
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA, all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic", analysis_type = "survival", predictor = predictor_name,
          num_cuts = num_cuts, criterion = criterion, covariates = covariates,
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
      if (!quiet) cli::cli_inform("Not enough data ({nrow(userdata)}) for nmin ({nmin}) and {num_cuts} cut(s). Returning NA cuts.")
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA, all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic", analysis_type = "survival", predictor = predictor_name,
          num_cuts = num_cuts, criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }
    stats_per_cut <- sapply(search_grid, .get_stat, num_cuts = 1, data_in = userdata,
                            criterion = criterion, cov_formula = cov_part, nmin = nmin,
                            fit_null = fit_null_model) # Pass null model

    if (all(is.na(stats_per_cut))) {
      if (!quiet) cli::cli_inform("No valid cut-points found due to model failures or constraints.")
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA, all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic", analysis_type = "survival", predictor = predictor_name,
          num_cuts = num_cuts, criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }
    best_idx <- if (direction == "min") which.min(stats_per_cut) else which.max(stats_per_cut)
    best_cut_val <- search_grid[best_idx]
    best_stat <- stats_per_cut[best_idx]
    all_stats_df <- data.frame(cut1 = search_grid, stat = stats_per_cut)

  } else { # num_cuts == 2
    if (!quiet) cli::cli_alert_info("Searching for 2 cuts is computationally intensive...")

    # Set up sequential backend
    foreach::registerDoSEQ()

    possible_c1_indices <- nmin:(nrow(userdata) - (2 * nmin))
    grid1_values <- unique(userdata$factor[possible_c1_indices])

    if (length(grid1_values) == 0) {
      if (!quiet) cli::cli_inform("Not enough data ({nrow(userdata)}) for nmin ({nmin}) and {num_cuts} cut(s). Returning NA cuts.")
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA, all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic", analysis_type = "survival", predictor = predictor_name,
          num_cuts = num_cuts, criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }

    results_list <- foreach::foreach(c1 = grid1_values, .combine = 'rbind', .export = c(".get_stat")) %do% {
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
        stat <- .get_stat(c(c1, c2), 2, userdata, criterion, cov_part, nmin, fit_null = fit_null_model) # Pass null model
        if (is.na(stat)) next

        is_better <- if (direction == "min") (stat < best_local_stat) else (stat > best_local_stat)
        if (is_better && !is.infinite(stat)) {
          best_local_stat <- stat
          best_local_c2 <- c2
        }
      }
      if (is.na(best_local_c2)) return(NULL)
      data.frame(stat = best_local_stat, c1 = c1, c2 = c2)
    }

    if (is.null(results_list) || nrow(results_list) == 0) {
      if (!quiet) cli::cli_inform("Not enough data ({nrow(userdata)}) for nmin ({nmin}) and {num_cuts} cut(s). Returning NA cuts.")
      # Return a valid na_result structure
      return(list(
        optimal_cuts = rep(NA, num_cuts), optimal_stat = NA, all_stats = NULL,
        userdata = userdata,
        parameters = list(
          method = "systematic", analysis_type = "survival", predictor = predictor_name,
          num_cuts = num_cuts, criterion = criterion, covariates = covariates,
          nmin = nmin, quiet = quiet
        )
      ))
    }

    best_idx <- if (direction == "min") which.min(results_list$stat) else which.max(results_list$stat)
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

# --- Internal Core Logic: Calculate Statistic (for systematic search) ---
.get_stat <- function(cuts, num_cuts, data_in, criterion, cov_formula, nmin, fit_null = NULL) {
  breaks <- sort(unique(c(-Inf, cuts, Inf)))
  num_intervals <- length(breaks) - 1
  data_in$group <- factor(cut(data_in$factor, breaks = breaks, labels = 1:num_intervals))

  if (any(table(data_in$group) < nmin) || nlevels(data_in$group) != (num_cuts + 1)) {
    return(NA)
  }

  formula_str <- paste("Surv(time, event) ~ group", cov_formula)

  if (criterion == "logrank") {
    # survdiff does not support covariates, so we must use a Cox model's score test.
    if (cov_formula == "") {
      fit <- tryCatch(survival::survdiff(as.formula(formula_str), data = data_in), error = function(e) NULL)
      if (is.null(fit)) return(NA)
      return(fit$chisq)
    } else {
      # Use the score test from the Cox model (equivalent to log-rank with covariates)
      fit <- tryCatch(survival::coxph(as.formula(formula_str), data = data_in), error = function(e) NULL)
      # *** CORRECTED: Use fit$score, not fit$logtest ***
      if (is.null(fit) || is.null(fit$score)) return(NA)
      return(fit$score)
    }
  } else { # Cox-based criteria
    fit <- tryCatch(survival::coxph(as.formula(formula_str), data = data_in), error = function(e) NULL)
    if (is.null(fit)) return(NA)

    if (criterion == "hazard_ratio") {
      if (is.null(fit$coefficients)) return(NA)
      # Optimized: get coef directly, don't use summary()
      # Note: group1 is ref, so we get coef for group2 (for num_cuts = 1)
      coef_name <- paste0("group", num_cuts + 1) # Generalize for potential future use
      if (!coef_name %in% names(fit$coefficients)) return(NA) # Ensure coef exists
      return(exp(fit$coefficients[coef_name]))
    } else if (criterion == "p_value") {
      # Optimized: Use Likelihood Ratio Test
      if (is.null(fit_null) || is.null(fit$loglik)) return(NA)
      loglik0 <- fit_null$loglik[2] # Loglik from null model
      loglik1 <- fit$loglik[2]      # Loglik from full model (with group)
      lrt_stat <- 2 * (loglik1 - loglik0)
      df_diff <- num_cuts # Df is number of groups - 1
      p_val <- stats::pchisq(lrt_stat, df = df_diff, lower.tail = FALSE)
      return(p_val)
    }
  }
}

# --- S3 Methods for Print, Summary, Plot ---

#' @rdname find_cutpoint
#' @export
print.find_cutpoint <- function(x, ...) {
  if (is.null(x) || any(is.na(x$optimal_cuts))) {
    cli::cli_inform("No optimal cut-point could be determined with the given parameters.")
    return(invisible(x))
  }

  method_name <- tools::toTitleCase(x$parameters$method)
  cli::cli_h1("Optimal Cut-point Analysis for Survival Data ({method_name})")

  stat_label <- switch(x$parameters$criterion,
                       "logrank" = "Optimal Log-Rank Statistic",
                       "hazard_ratio" = "Optimal Hazard Ratio",
                       "p_value" = "Optimal P-value")

  stat_val <- x$optimal_stat
  stat_val_fmt <- round(stat_val, 4)

  # For p_value (genetic), the stat is the LRT statistic, not the p-value itself.
  # Clarify this in the print output.
  if (x$parameters$criterion == "p_value" && x$parameters$method == "genetic") {
    stat_label <- "Optimal LRT Statistic"
  }


  cli::cli_bullets(c(
    "*" = "Predictor: {.strong {x$parameters$predictor}}",
    "*" = "Criterion: {.strong {x$parameters$criterion}}",
    "*" = "{stat_label}: {.strong {stat_val_fmt}}",
    "v" = "Recommended Cut-point(s): {.strong {paste(round(x$optimal_cuts, 3), collapse = ', ')}}"
  ))
  invisible(x)
}

#' @rdname find_cutpoint
#' @export
summary.find_cutpoint <- function(object, show_model = TRUE, show_group_counts = TRUE,
                                  show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE, ...) {
  cli::cli_h1("Optimal Cut-point Analysis for Survival Data ({tools::toTitleCase(object$parameters$method)})")

  if (is.null(object) || any(is.na(object$optimal_cuts))) {
    cli::cli_inform("No valid optimal cut-point was found with the given parameters.")
    return(invisible(object))
  }

  data <- object$userdata
  cuts <- object$optimal_cuts
  num_cuts <- object$parameters$num_cuts

  data$group <- factor(cut(data$factor, breaks = c(-Inf, cuts, Inf),
                           labels = paste0("G", 1:(num_cuts + 1))))

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

  formula_str <- "survival::Surv(time, event) ~ group"
  if (!is.null(object$parameters$covariates)) {
    formula_str <- paste(formula_str, "+", paste(object$parameters$covariates, collapse = " + "))
  }
  fit_cox <- tryCatch(survival::coxph(as.formula(formula_str), data = data), error = function(e) NULL)

  if (show_model) {
    cli::cli_h2("Final Cox Model Summary")
    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model for summary: model convergence failed.")
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
      param_bullets <- c(param_bullets, "*" = "Covariates: {paste(params$covariates, collapse = ', ')}")
    }
    cli::cli_bullets(param_bullets)
  }

  invisible(object)
}

#' @rdname find_cutpoint
#' @export
plot.find_cutpoint <- function(x, type = "outcome", reference_group = NULL, ...) {
  type <- match.arg(type, choices = c("outcome", "distribution", "forest"))

  if (is.null(x) || any(is.na(x$optimal_cuts))) {
    cli::cli_inform("Cannot generate plot, no valid cut-point found.")
    return(invisible(NULL))
  }

  cuts <- x$optimal_cuts
  data <- x$userdata
  num_cuts <- length(cuts)
  data$group <- factor(cut(data$factor, breaks = c(-Inf, cuts, Inf),
                           labels = paste0("G", 1:(num_cuts + 1))))

  if (type == "distribution") {
    p <- ggplot2::ggplot(x$userdata, ggplot2::aes(x = .data$factor)) +
      ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)), bins = 30, fill = "#56B4E9", color = "black", alpha = 0.7) +
      ggplot2::geom_density(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_vline(xintercept = cuts, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
      ggplot2::labs(title = "Distribution of Predictor with Optimal Cut-points", x = x$parameters$predictor, y = "Density") +
      ggplot2::theme_minimal()
    return(p)

  } else if (type == "outcome") {
    fit <- survival::survfit(Surv(time, event) ~ group, data = data)
    p <- survminer::ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE,
                               legend.title = "Groups",
                               palette = "jco",
                               ggtheme = ggplot2::theme_minimal())
    p$plot <- p$plot + ggplot2::labs(title = paste("Survival Curves by", x$parameters$predictor, "Group"))
    return(p)

  } else if (type == "forest") {
    if (!requireNamespace("broom", quietly = TRUE)) {
      stop("Package 'broom' is required for the forest plot. Please install it.", call. = FALSE)
    }

    cli::cli_alert_info("Generating Forest Plot of Hazard Ratios...")

    group_levels <- levels(data$group)
    if (is.null(reference_group) || !(reference_group %in% group_levels)) {
      reference_group <- group_levels[1]
      cli::cli_alert_warning("Reference group not specified or invalid, defaulting to: {reference_group}")
    }
    data$group <- stats::relevel(data$group, ref = reference_group)

    formula_str <- "survival::Surv(time, event) ~ group"
    if (!is.null(x$parameters$covariates)) {
      formula_str <- paste(formula_str, "+", paste(x$parameters$covariates, collapse = " + "))
    }
    fit_cox <- tryCatch(survival::coxph(as.formula(formula_str), data = data), error = function(e) NULL)

    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model for forest plot: model convergence failed.")
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

