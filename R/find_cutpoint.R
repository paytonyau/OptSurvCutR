#' Find Optimal Cut-points with Advanced Features
#'
#' @description
#' Determines optimal cut-point(s) using a systematic search (for 1-2 cuts)
#' or a genetic algorithm (for any number of cuts). The genetic algorithm now
#' supports both survival and logistic regression outcomes.
#'
#' @param data A data frame.
#' @param predictor The name of the predictor variable.
#' @param num_cuts The number of cut-points to find.
#' @param method The algorithm to use: "systematic" or "genetic".
#' @param covariates Covariates to adjust for.
#' @param outcome_time Time variable for survival analysis.
#' @param outcome_event Event variable for survival analysis.
#' @param outcome_binary Binary outcome variable.
#' @param nmin Minimum number of observations per group.
#' @param ... Additional arguments for the genetic algorithm (e.g., `popSize`, `maxiter`).
#'
#' @return An object containing the results of the analysis.
#' @importFrom stats na.omit as.formula pchisq glm binomial predict confint dnorm uniroot quantile sd fitted AIC rnorm runif
#' @importFrom survival coxph Surv survdiff
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success cli_bullets cli_progress_bar cli_progress_update cli_progress_done
#' @importFrom pROC roc auc coords ggroc
#' @importFrom ggplot2 ggplot aes .data geom_line geom_vline labs theme_minimal geom_histogram geom_density
#' @importFrom foreach %dopar%
#' @export
find_cutpoint <- function(data, predictor, num_cuts = 1, method = "systematic", covariates = NULL,
                          outcome_time = NULL, outcome_event = NULL,
                          outcome_binary = NULL, nmin = 20, ...) {

  # --- 1. Input Validation and Data Prep ---
  method <- match.arg(method, choices = c("systematic", "genetic"))
  if (is.null(predictor)) stop("A 'predictor' variable must be specified.", call. = FALSE)

  if (!is.null(outcome_time) && !is.null(outcome_event)) {
    analysis_type <- "survival"
  } else if (!is.null(outcome_binary)) {
    analysis_type <- "logistic"
  } else {
    stop("You must specify outcome variables.", call. = FALSE)
  }

  required_vars <- c(predictor, covariates, outcome_time, outcome_event, outcome_binary)
  userdata <- data[, unique(required_vars), drop = FALSE]
  userdata <- stats::na.omit(userdata)

  names(userdata)[names(userdata) == predictor] <- "factor"
  if(analysis_type == "survival") {
    names(userdata)[names(userdata) == outcome_time] <- "time"
    names(userdata)[names(userdata) == outcome_event] <- "event"
  } else {
    names(userdata)[names(userdata) == outcome_binary] <- "outcome"
  }

  # --- 2. Route to appropriate method ---
  if (method == "systematic") {
    # Bundle parameters into a list to avoid "unused argument" warnings
    params <- list(userdata = userdata, num_cuts = num_cuts,
                   analysis_type = analysis_type, covariates = covariates, nmin = nmin)
    .systematic_search(params)
  } else { # genetic
    .genetic_search(userdata = userdata, num_cuts = num_cuts,
                    analysis_type = analysis_type, covariates = covariates, nmin = nmin, ...)
  }
}

# --- Internal Helper: Systematic Search ---
.systematic_search <- function(p) { # p is for parameters
  if (!p$num_cuts %in% c(1, 2)) {
    stop("Systematic search only supports num_cuts = 1 or 2.", call. = FALSE)
  }

  cli::cli_alert_info("Starting systematic search for {p$num_cuts} optimal cut-point(s)...")

  p$userdata <- p$userdata[order(p$userdata$factor), ]
  cov_part <- if (!is.null(p$covariates)) paste(" +", paste(p$covariates, collapse = " + ")) else ""

  best_stat <- if(p$analysis_type == "survival") Inf else -Inf
  best_cut_val <- rep(NA, p$num_cuts)
  all_stats_df <- NULL

  if (p$num_cuts == 1) {
    search_grid <- unique(p$userdata$factor[p$nmin:(nrow(p$userdata) - p$nmin)])
    stats_per_cut <- sapply(search_grid, .get_stat, num_cuts = 1, data_in = p$userdata, analysis = p$analysis_type, cov_formula = cov_part, nmin = p$nmin)
    if(all(is.na(stats_per_cut))) stop("Could not find any valid cut-points.", call. = FALSE)

    best_idx <- if(p$analysis_type == "survival") which.min(stats_per_cut) else which.max(stats_per_cut)
    best_cut_val <- search_grid[best_idx]
    best_stat <- stats_per_cut[best_idx]
    all_stats_df <- data.frame(cut1 = search_grid, stat = stats_per_cut)

  } else { # num_cuts == 2
    cli::cli_alert_info("Searching for 2 cuts is computationally intensive and may be slow.")

    grid1 <- unique(p$userdata$factor[p$nmin:(nrow(p$userdata) - 2 * p$nmin)])
    pb <- cli::cli_progress_bar("Evaluating cut-point pairs", total = length(grid1))

    for (c1 in grid1) {
      cli::cli_progress_update(id = pb)
      min_idx_c2 <- which(p$userdata$factor > c1)[p$nmin]
      if (is.na(min_idx_c2)) next
      max_idx_c2 <- nrow(p$userdata) - p$nmin
      if (min_idx_c2 > max_idx_c2) next

      grid2 <- unique(p$userdata$factor[min_idx_c2:max_idx_c2])

      for (c2 in grid2) {
        stat <- .get_stat(c(c1, c2), num_cuts = 2, data_in = p$userdata, analysis = p$analysis_type, cov_formula = cov_part, nmin = p$nmin)
        if (is.na(stat)) next

        is_better <- if(p$analysis_type == "survival") (stat < best_stat) else (stat > best_stat)
        if (is_better) {
          best_stat <- stat
          best_cut_val <- c(c1, c2)
        }
      }
    }
    cli::cli_progress_done(id = pb)
  }

  output <- list(
    best_cut = best_cut_val,
    best_stat = best_stat,
    all_stats = all_stats_df,
    parameters = list(analysis_type = p$analysis_type, predictor = "factor", num_cuts = p$num_cuts, covariates = p$covariates, nmin = p$nmin),
    userdata = p$userdata
  )
  class(output) <- "find_cutpoint_systematic"
  cli::cli_alert_success("Systematic search complete.")
  return(output)
}

# --- Internal Helper: Genetic Algorithm ---
.genetic_search <- function(userdata, num_cuts, analysis_type, covariates, nmin, ...) {
  cli::cli_alert_info("Starting genetic algorithm for {num_cuts} cut-point(s)...")

  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

  if (analysis_type == "survival") {
    fitness_function <- function(cuts, data) {
      val <- .get_stat(cuts, num_cuts, data, "survival", cov_part, nmin)
      if(is.na(val)) return(Inf)
      return(val)
    }
    direction <- "min"
  } else { # Logistic
    fitness_function <- function(cuts, data) {
      val <- .get_stat(cuts, num_cuts, data, "logistic", cov_part, nmin)
      if(is.na(val)) return(Inf)
      return(val)
    }
    direction <- "min" # We minimize AIC
  }

  ga_result <- .run_ga(fitness_function, direction = direction, num_cuts = num_cuts, data = userdata, ...)

  output <- list(
    optimal_cuts = sort(ga_result$best_solution),
    userdata = userdata,
    parameters = list(analysis_type = analysis_type, num_cuts = num_cuts, covariates = covariates, nmin = nmin)
  )
  class(output) <- "find_cutpoint_genetic"
  cli::cli_alert_success("Genetic algorithm complete.")
  return(output)
}

#' Perform Permutation Testing for a Systematic Cut-point Search
#'
#' @description
#' Calculates an adjusted p-value for a result from a systematic `find_cutpoint`
#' search to correct for multiple testing.
#'
#' @param cutpoint_result An object from `find_cutpoint(method = "systematic")`.
#' @param permutations The number of permutations to perform.
#'
#' @return A list containing the original statistic and the adjusted p-value.
#' @export
permute_cutpoint <- function(cutpoint_result, permutations = 1000) {
  if (!inherits(cutpoint_result, "find_cutpoint_systematic")) {
    stop("Permutation testing is only applicable to results from a systematic search.", call. = FALSE)
  }

  params <- cutpoint_result$parameters
  observed_stat <- cutpoint_result$best_stat

  cli::cli_alert_info("Running {permutations} permutations for p-value adjustment...")
  pb <- cli::cli_progress_bar("Permutations", total = permutations)

  perm_stats <- replicate(permutations, {
    cli::cli_progress_update(id = pb)
    perm_data <- cutpoint_result$userdata
    perm_data$factor <- sample(perm_data$factor)

    # Bundle params for the internal call
    p_perm <- list(userdata = perm_data, num_cuts = params$num_cuts,
                   analysis_type = params$analysis_type,
                   covariates = params$covariates, nmin = params$nmin)

    perm_res <- .systematic_search(p_perm)
    return(perm_res$best_stat)
  })
  cli::cli_progress_done(id = pb)

  if(params$analysis_type == "survival") {
    p_adj <- mean(perm_stats <= observed_stat, na.rm = TRUE)
  } else {
    p_adj <- mean(perm_stats >= observed_stat, na.rm = TRUE)
  }

  cli::cli_alert_success("Permutation test complete.")
  return(list(original_statistic = observed_stat, adjusted_p_value = p_adj))
}


# --- S3 Methods for Print, Summary, Plot ---
#' @param x An object from `find_cutpoint`.
#' @param ... Additional arguments.
#' @rdname find_cutpoint
#' @export
print.find_cutpoint_systematic <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Analysis (Systematic)")
  cli::cli_bullets(c(
    "\u2714" = "Final Recommended Cut-point(s): {.strong {paste(round(x$best_cut, 3), collapse = ', ')}}"
  ))
  invisible(x)
}

#' @rdname find_cutpoint
#' @export
print.find_cutpoint_genetic <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Analysis (Genetic)")
  cli::cli_bullets(c(
    "\u2714" = "Final Recommended Cut-point(s): {.strong {paste(round(x$optimal_cuts, 3), collapse = ', ')}}"
  ))
  invisible(x)
}

#' @rdname find_cutpoint
#' @param object An object from `find_cutpoint`.
#' @export
summary.find_cutpoint_systematic <- function(object, ...) {
  cli::cli_h1("Summary for Systematic Search")
  print(object$parameters)
}

#' @rdname find_cutpoint
#' @export
summary.find_cutpoint_genetic <- function(object, ...) {
  cli::cli_h1("Summary for Genetic Algorithm Search")
  print(object$parameters)
}

#' @rdname find_cutpoint
#' @param type The type of plot.
#' @importFrom survminer ggsurvplot
#' @export
plot.find_cutpoint_systematic <- function(x, type = "outcome", ...) {
  plot.find_cutpoint_genetic(x, type = type, ...)
}

#' @rdname find_cutpoint
#' @export
plot.find_cutpoint_genetic <- function(x, type = "outcome", ...) {
  type <- match.arg(type, choices = c("outcome", "distribution"))
  cuts <- if(!is.null(x$best_cut)) x$best_cut else x$optimal_cuts

  if (type == "distribution") {
    p <- ggplot2::ggplot(x$userdata, ggplot2::aes(x = .data$factor)) +
      ggplot2::geom_histogram(aes(y = ggplot2::after_stat(density)), bins = 30, fill = "grey", color = "white") +
      ggplot2::geom_vline(xintercept = cuts, color = "red", linetype = "dashed", linewidth = 1.2) +
      ggplot2::labs(title = "Distribution of Predictor with Optimal Cut-points")
    return(p)
  }

  if (type == "outcome") {
    data <- x$userdata
    data$group <- factor(cut(data$factor, breaks = c(-Inf, cuts, Inf)))
    if (x$parameters$analysis_type == "survival") {
      fit <- survival::survfit(Surv(time, event) ~ group, data = data)
      p <- survminer::ggsurvplot(fit, data = data, pval = TRUE, risk.table = TRUE)
      return(p)
    } else {
      roc_obj <- pROC::roc(data$outcome, data$factor, quiet = TRUE)
      p <- pROC::ggroc(roc_obj) + ggplot2::labs(title = "ROC Curve for Predictor")
      return(p)
    }
  }
}

# --- Internal Core Logic ---
.get_stat <- function(cuts, num_cuts, data_in, analysis, cov_formula, nmin) {
  breaks <- c(-Inf, sort(cuts), Inf)
  data_in$group <- factor(cut(data_in$factor, breaks = breaks))

  if(any(table(data_in$group) < nmin) || nlevels(data_in$group) < (num_cuts + 1)) return(NA)

  formula_str <- paste("~ group", cov_formula)

  if(analysis == "survival") {
    formula_str <- paste("Surv(time, event)", formula_str)
    fit <- tryCatch(survival::coxph(as.formula(formula_str), data = data_in), error = function(e) NULL)
    if(is.null(fit)) return(NA)
    return(summary(fit)$sctest["pvalue"])
  } else { # logistic
    formula_str <- paste("outcome", formula_str)
    fit <- tryCatch(stats::glm(as.formula(formula_str), data = data_in, family = "binomial"), error=function(e) NULL)
    if(is.null(fit)) return(NA)
    return(stats::AIC(fit))
  }
}

.run_ga <- function(fitness, direction, num_cuts, data, popSize = 50, maxiter = 100, ...) {
  min_val <- min(data$factor)
  max_val <- max(data$factor)

  population <- t(replicate(popSize, sort(stats::runif(num_cuts, min_val, max_val))))

  best_fitness <- if(direction == "min") Inf else -Inf
  best_solution <- population[1, ]

  for (i in 1:maxiter) {
    fitness_scores <- apply(population, 1, fitness, data = data)

    current_best_idx <- if(direction == "min") which.min(fitness_scores) else which.max(fitness_scores)
    current_best_fitness <- fitness_scores[current_best_idx]

    if ((direction == "min" && current_best_fitness < best_fitness) ||
        (direction == "max" && current_best_fitness > best_fitness)) {
      best_fitness <- current_best_fitness
      best_solution <- population[current_best_idx, ]
    }

    parents_idx <- sample(1:popSize, size = popSize, replace = TRUE)
    parents <- population[parents_idx, ]

    for (j in 1:popSize) {
      if (stats::runif(1) < 0.8 && num_cuts > 1) {
        parent2_idx <- sample(1:popSize, 1)
        crossover_point <- sample(1:(num_cuts - 1), 1)
        population[j, ] <- c(parents[j, 1:crossover_point], parents[parent2_idx, (crossover_point+1):num_cuts])
      }
      if (stats::runif(1) < 0.1) {
        mutate_point <- sample(1:num_cuts, 1)
        population[j, mutate_point] <- population[j, mutate_point] + stats::rnorm(1, 0, sd(data$factor)/10)
      }
    }
  }
  return(list(best_solution = best_solution, best_fitness = best_fitness))
}

# Suppress NOTE about ggplot2 'density' aesthetic
utils::globalVariables(c("density"))
