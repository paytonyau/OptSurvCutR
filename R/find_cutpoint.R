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
#' @param seed An optional integer to set the random seed for reproducible results
#'   when `method = "genetic"`.
#' @param maxiter The number of generations for the genetic algorithm. Default is 100.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores for the systematic search.
#' @param quiet Logical. If TRUE, suppresses the printing of the final result object.
#' @param ... Additional arguments passed to the genetic algorithm (e.g., `popSize`).
#' @param x An object from `find_cutpoint`.
#' @param object An object from `find_cutpoint`.
#' @param show_model Logical. If TRUE, shows the full summary of the final Cox model.
#' @param show_group_counts Logical. If TRUE, shows the number of subjects and events in each group.
#' @param show_medians Logical. If TRUE, shows the median survival for each group.
#' @param show_ph_test Logical. If TRUE, shows the proportional hazards assumption test.
#' @param show_params Logical. If TRUE, shows the parameters of the original function call.
#' @param type The type of plot to generate: "outcome" (a survival plot) or "distribution".
#'
#' @return An object of class `find_cutpoint` containing the optimal cut-points,
#'   the corresponding statistic, and other parameters used in the analysis.
#' @importFrom stats na.omit as.formula pchisq sd rnorm runif anova aggregate
#' @importFrom survival Surv survfit survdiff coxph cox.zph
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_bullets cli_progress_bar cli_progress_update cli_progress_done cli_alert_danger cli_warn
#' @importFrom ggplot2 ggplot aes .data geom_line geom_vline labs theme_minimal geom_histogram geom_density
#' @importFrom foreach %dopar%
#' @importFrom doParallel registerDoParallel
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom survminer ggsurvplot
#' @importFrom tools toTitleCase
#' @export
find_cutpoint <- function(data, predictor, outcome_time, outcome_event, num_cuts = 1, method = "systematic",
                          criterion = "logrank", covariates = NULL, nmin = 20, seed = NULL, maxiter = 100, use_parallel = TRUE, quiet = FALSE, ...) {

  # --- 1. Input Validation and Data Prep ---
  method <- match.arg(method, choices = c("systematic", "genetic"))
  criterion <- match.arg(criterion, choices = c("logrank", "hazard_ratio", "p_value"))

  # --- FIXED: Consolidated and robust input validation for all required columns ---
  if (is.null(predictor)) stop("A 'predictor' variable must be specified.", call. = FALSE)
  if (is.null(outcome_time) || is.null(outcome_event)) stop("Both 'outcome_time' and 'outcome_event' must be specified.", call. = FALSE)

  required_vars <- c(predictor, outcome_time, outcome_event, covariates)
  missing_vars <- setdiff(required_vars, names(data))

  if (length(missing_vars) > 0) {
    stop(paste0("The following required column(s) were not found in the data: '",
                paste(missing_vars, collapse = "', '"), "'"), call. = FALSE)
  }
  # --- End of fix ---

  if (criterion == "hazard_ratio" && num_cuts > 1) {
    stop("'hazard_ratio' criterion is only supported for num_cuts = 1.", call. = FALSE)
  }

  if (method == "genetic" && !is.null(seed)) {
    set.seed(seed)
  }

  userdata <- data[, unique(required_vars), drop = FALSE]
  userdata <- stats::na.omit(userdata)

  original_predictor_name <- predictor
  names(userdata)[names(userdata) == predictor] <- "factor"
  names(userdata)[names(userdata) == outcome_time] <- "time"
  names(userdata)[names(userdata) == outcome_event] <- "event"

  if (any(userdata$time < 0)) {
    stop("Time variable must be non-negative.", call. = FALSE)
  }

  if (nmin > 0 && nmin < 1) nmin <- floor(nmin * nrow(userdata))
  if (nrow(userdata) < nmin * (num_cuts + 1)) {
    stop(paste0("Not enough data (", nrow(userdata), ") for nmin (", nmin, ") and ", num_cuts, " cut(s)."), call. = FALSE)
  }

  # --- 2. Route to appropriate search method ---
  params <- list(userdata = userdata, num_cuts = num_cuts, criterion = criterion,
                 covariates = covariates, nmin = nmin, predictor_name = original_predictor_name,
                 maxiter = maxiter, use_parallel = use_parallel, ...)

  output <- if (method == "systematic") {
    do.call(.systematic_search, params)
  } else { # genetic
    cli::cli_alert_info("Starting genetic algorithm for {num_cuts} cut-point(s) using '{criterion}' criterion...")

    ga_result <- .run_genetic_search(
      target = userdata$factor,
      numcut = num_cuts,
      time = userdata$time,
      censor = userdata$event,
      confound = if(!is.null(covariates)) userdata[, covariates, drop=FALSE] else NULL,
      nmin = nmin,
      criterion = criterion,
      numgen = maxiter,
      ...
    )

    if (is.null(ga_result) || !is.finite(ga_result$value)) {
      cli::cli_warn("Genetic algorithm could not find a valid solution with the given constraints.")
      optimal_cuts <- rep(NA, num_cuts)
      optimal_stat <- NA
    } else {
      optimal_cuts <- sort(ga_result$par)
      optimal_stat <- ga_result$value
      if (criterion == "p_value") {
        optimal_stat <- 1 - optimal_stat
      }
    }

    list(
      optimal_cuts = optimal_cuts,
      optimal_stat = optimal_stat,
      all_stats = NULL,
      userdata = userdata,
      parameters = list(method = "genetic", analysis_type = "survival", predictor = original_predictor_name,
                        num_cuts = num_cuts, criterion = criterion,
                        covariates = covariates, nmin = nmin)
    )
  }

  class(output) <- "find_cutpoint"

  # Only print the output if not in quiet mode
  if (!quiet) {
    print(output)
  }

  invisible(output)
}


# --- Internal Helper: Systematic Search ---
.systematic_search <- function(userdata, num_cuts, criterion, covariates, nmin, predictor_name, use_parallel, ...) {
  if (!num_cuts %in% c(1, 2)) {
    stop("Systematic search currently only supports num_cuts = 1 or 2.", call. = FALSE)
  }
  cli::cli_alert_info("Starting systematic search for {num_cuts} optimal cut-point(s) using '{criterion}' criterion...")
  userdata <- userdata[order(userdata$factor), ]
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

  direction <- if (criterion == "p_value") "min" else "max"
  best_stat <- if (direction == "min") Inf else -Inf
  best_cut_val <- rep(NA, num_cuts)
  all_stats_df <- NULL

  if (num_cuts == 1) {
    search_grid <- unique(userdata$factor[nmin:(nrow(userdata) - nmin)])
    if(length(search_grid) == 0) {
      best_stat <- NA # Indicate failure
    } else {
      stats_per_cut <- sapply(search_grid, .get_stat, num_cuts = 1, data_in = userdata,
                              criterion = criterion, cov_formula = cov_part, nmin = nmin)
      if (all(is.na(stats_per_cut))) {
        best_stat <- NA
      } else {
        best_idx <- if (direction == "min") which.min(stats_per_cut) else which.max(stats_per_cut)
        best_cut_val <- search_grid[best_idx]
        best_stat <- stats_per_cut[best_idx]
        all_stats_df <- data.frame(cut1 = search_grid, stat = stats_per_cut)
      }
    }
  } else { # num_cuts == 2
    cli::cli_alert_info("Searching for 2 cuts is computationally intensive. Using parallel processing if available.")

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

    possible_c1_indices <- nmin:(nrow(userdata) - (2 * nmin))
    grid1_values <- unique(userdata$factor[possible_c1_indices])

    results_list <- foreach::foreach(c1 = grid1_values, .combine = 'rbind', .export = c(".get_stat")) %dopar% {
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

      for (c2 in grid2_values) {
        stat <- .get_stat(c(c1, c2), 2, userdata, criterion, cov_part, nmin)
        if (is.na(stat)) next

        is_better <- if (direction == "min") (stat < best_local_stat) else (stat > best_local_stat)
        if (is_better && !is.infinite(stat)) {
          best_local_stat <- stat
          best_local_c2 <- c2
        }
      }
      if(is.na(best_local_c2)) return(NULL)
      data.frame(stat = best_local_stat, c1 = c1, c2 = best_local_c2)
    }

    if(!is.null(results_list) && nrow(results_list) > 0){
      best_idx <- if(direction == "min") which.min(results_list$stat) else which.max(results_list$stat)
      best_stat <- results_list$stat[best_idx]
      best_cut_val <- c(results_list$c1[best_idx], results_list$c2[best_idx])
    } else {
      best_stat <- NA
    }
  }

  if(is.na(best_stat) || is.infinite(best_stat)){
    cli::cli_warn("Systematic search could not find any valid cut-points with the given constraints.")
    return(list(optimal_cuts = rep(NA, num_cuts), optimal_stat = NA, all_stats = NULL, userdata = userdata,
                parameters = list(method = "systematic", analysis_type = "survival", predictor = predictor_name,
                                  num_cuts = num_cuts, criterion = criterion,
                                  covariates = covariates, nmin = nmin)))
  }

  output <- list(
    optimal_cuts = best_cut_val,
    optimal_stat = best_stat,
    all_stats = all_stats_df,
    parameters = list(method = "systematic", analysis_type = "survival", predictor = predictor_name,
                      num_cuts = num_cuts, criterion = criterion,
                      covariates = covariates, nmin = nmin),
    userdata = userdata
  )
  cli::cli_alert_success("Systematic search complete.")
  return(output)
}


# --- Internal Core Logic: Calculate Statistic (for systematic search) ---
.get_stat <- function(cuts, num_cuts, data_in, criterion, cov_formula, nmin) {
  breaks <- sort(unique(c(-Inf, cuts, Inf)))
  num_intervals <- length(breaks) - 1
  data_in$group <- factor(cut(data_in$factor, breaks = breaks, labels = 1:num_intervals))

  if (any(table(data_in$group) < nmin) || nlevels(data_in$group) != (num_cuts + 1)) {
    return(NA)
  }

  formula_str <- paste("Surv(time, event) ~ group", cov_formula)

  if (criterion == "logrank") {
    fit <- tryCatch(survival::survdiff(as.formula(formula_str), data = data_in), error = function(e) NULL)
    if (is.null(fit)) return(NA)
    return(fit$chisq)
  } else { # Cox-based criteria
    fit <- tryCatch(survival::coxph(as.formula(formula_str), data = data_in), error = function(e) NULL)
    if (is.null(fit)) return(NA)

    if (criterion == "hazard_ratio") {
      return(summary(fit)$conf.int[1, "exp(coef)"])
    } else if (criterion == "p_value") {
      return(summary(fit)$logtest["pvalue"])
    }
  }
}

# --- S3 Methods for Print, Summary, Plot ---

#' @rdname find_cutpoint
#' @export
print.find_cutpoint <- function(x, ...) {
  if(is.null(x) || any(is.na(x$optimal_cuts))) {
    cli::cli_alert_danger("No optimal cut-point could be determined with the given parameters.")
    return(invisible(x))
  }

  method_name <- tools::toTitleCase(x$parameters$method)
  cli::cli_h1("Optimal Cut-point Analysis for Survival Data ({method_name})")

  stat_label <- switch(x$parameters$criterion,
                       "logrank" = "Optimal Log-Rank Statistic",
                       "hazard_ratio" = "Optimal Hazard Ratio",
                       "p_value" = "Optimal P-value")

  cli::cli_bullets(c(
    "*" = "Predictor: {.strong {x$parameters$predictor}}",
    "*" = "Criterion: {.strong {x$parameters$criterion}}",
    "*" = "{stat_label}: {.strong {round(x$optimal_stat, 4)}}",
    "v" = "Recommended Cut-point(s): {.strong {paste(round(x$optimal_cuts, 3), collapse = ', ')}}"
  ))
  invisible(x)
}

#' @rdname find_cutpoint
#' @export
summary.find_cutpoint <- function(object, show_model = TRUE, show_group_counts = TRUE,
                                  show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE, ...) {

  if(is.null(object) || any(is.na(object$optimal_cuts))) {
    cli::cli_alert_danger("Cannot generate summary because no optimal cut-point was found.")
    return(invisible(NULL))
  }

  print(object)

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
  fit_cox <- survival::coxph(as.formula(formula_str), data = data)

  if (show_model) {
    cli::cli_h2("Final Cox Model Summary")
    print(summary(fit_cox))
  }

  if (show_ph_test) {
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
    if(!is.null(params$covariates)){
      param_bullets <- c(param_bullets, "*" = "Covariates: {paste(params$covariates, collapse=', ')}")
    }
    cli::cli_bullets(param_bullets)
  }

  invisible(object)
}


#' @rdname find_cutpoint
#' @export
plot.find_cutpoint <- function(x, type = "outcome", ...) {
  type <- match.arg(type, choices = c("outcome", "distribution"))

  if(is.null(x) || any(is.na(x$optimal_cuts))) {
    cli::cli_alert_warning("Cannot generate plot, no valid cut-point found.")
    return(invisible(NULL))
  }

  cuts <- x$optimal_cuts

  if (type == "distribution") {
    p <- ggplot2::ggplot(x$userdata, ggplot2::aes(x = .data$factor)) +
      ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)), bins = 30, fill = "#56B4E9", color = "black", alpha = 0.7) +
      ggplot2::geom_density(color = "#0072B2", linewidth = 1) +
      ggplot2::geom_vline(xintercept = cuts, color = "#D55E00", linetype = "dashed", linewidth = 1.2) +
      ggplot2::labs(title = "Distribution of Predictor with Optimal Cut-points", x = x$parameters$predictor, y = "Density") +
      ggplot2::theme_minimal()
    return(p)
  }

  if (type == "outcome") {
    data <- x$userdata
    data$group <- factor(cut(data$factor, breaks = c(-Inf, cuts, Inf)))
    fit <- survival::survfit(Surv(time, event) ~ group, data = data)
    p <- survminer::ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE,
                               legend.title = "Groups",
                               palette = "jco",
                               ggtheme = ggplot2::theme_minimal())
    p$plot <- p$plot + ggplot2::labs(title = paste("Survival Curves by", x$parameters$predictor, "Group"))
    return(p)
  }
}
