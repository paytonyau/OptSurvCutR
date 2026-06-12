# ===================================================================
# INTERNAL ENGINE: GENETIC SEARCH (rgenoud)
# Handles evolutionary threshold discovery for multi-dimensional cuts
# ===================================================================

#' Internal helper: Objective Function for Genetic Algorithm
#'
#' @description
#' Called by `rgenoud::genoud`. Translates incoming discrete integer board indices
#' back to regularized biomarker coordinates and calculates fitness (maximised).
#'
#' @inheritParams find_cutpoint
#' @inheritParams find_cutpoint_number
#' @param params Numeric vector. First `numcut` are integer indices of the grid; remaining are betas.
#' @param time Survival time vector.
#' @param censor Survival event vector.
#' @param target Continuous predictor vector.
#' @param confound Optional covariate data frame.
#' @param numcut Number of cut-points.
#' @param gap Minimum distance between cut-points (applied to decoded values).
#' @param loglik0 Log-likelihood of the null model.
#' @param cache Optional environment to cache and retrieve evaluations.
#' @param base_df Pre-allocated data.frame template for speed.
#' @param precompiled_formula Pre-parsed survival formula.
#' @param grid_pool The underlying regularized numeric coordinates vector.
#'
#' @return Single numeric fitness value.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.3} `na.rm = TRUE` used in quantile calculations.
#' @srrstats {G2.4a} `make.names()` sanitises covariate names.
#' @srrstats {G2.4d} Explicit conversion to factor used for group creation.
#' @srrstats {G2.13} Invalid gaps or insufficient N return `-Inf`.
#' @srrstats {G5.8} Edge cases (zero gap, infinite domain) return `-Inf`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival survdiff coxph
#' @importFrom stats pchisq
#' @noRd
.obj <- function(params, time, censor, target, confound, numcut, gap, nmin,
                 criterion, loglik0 = NA_real_, cache = NULL,
                 base_df = NULL, precompiled_formula = NULL, grid_pool = NULL) {
  # Decode incoming parameter indices back into numeric cuts
  winning_indices <- round(params[1:numcut])
  cutoff <- sort(grid_pool[winning_indices])

  if (!is.null(cache)) {
    key <- paste(round(cutoff, 6), collapse = "_")
    if (exists(key, envir = cache, inherits = FALSE)) {
      return(cache[[key]])
    }
  }

  if (length(unique(time)) <= 1) {
    return(-Inf)
  }

  if (numcut > 1 && min(diff(cutoff)) < gap) {
    return(-Inf)
  }
  if (cutoff[1] <= min(target, na.rm = TRUE) || cutoff[numcut] >= max(target, na.rm = TRUE)) {
    return(-Inf)
  }

  cut_design <- as.factor(findInterval(target, cutoff, left.open = TRUE) + 1L)
  if (length(table(cut_design)) != (numcut + 1) || min(table(cut_design)) < nmin) {
    return(-Inf)
  }

  data_for_fit <- base_df
  data_for_fit$cut_design <- cut_design

  stat_value <- switch(criterion,
    "logrank" = {
      if (is.null(confound) || ncol(confound) == 0) {
        fit <- tryCatch(survival::survdiff(precompiled_formula, data = data_for_fit), error = function(e) NULL)
        if (is.null(fit)) -Inf else fit$chisq
      } else {
        fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
        if (is.null(fit) || is.null(fit$score)) -Inf else fit$score
      }
    },
    "hazard_ratio" = {
      fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
      if (is.null(fit) || isTRUE(fit$nevent == 0)) {
        -Inf
      } else {
        target_coef <- paste0("cut_design", numcut + 1)
        if (!(target_coef %in% names(fit$coefficients)) || is.na(fit$coefficients[target_coef])) -Inf else exp(fit$coefficients[target_coef])
      }
    },
    "p_value" = {
      fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
      if (is.null(fit) || is.null(fit$loglik) || isTRUE(fit$nevent == 0)) {
        -Inf
      } else {
        if (!is.na(loglik0)) {
          pval <- stats::pchisq(-2 * (loglik0 - fit$loglik[2]), numcut, lower.tail = FALSE)
          1 - pval
        } else {
          sfit <- tryCatch(summary(fit), error = function(e) NULL)
          if (is.null(sfit) || is.null(sfit$logtest)) -Inf else 1 - sfit$logtest["pvalue"]
        }
      }
    },
    "loglik" = {
      beta_init <- params[(numcut + 1):length(params)]
      if (length(beta_init) == 0 || is.null(confound) || ncol(confound) == 0) {
        fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit), error = function(e) NULL)
        if (is.null(fit) || is.null(fit$loglik) || isTRUE(fit$nevent == 0)) -Inf else fit$loglik[2]
      } else {
        fit <- tryCatch(survival::coxph(precompiled_formula, data = data_for_fit, init = beta_init, iter.max = 0), error = function(e) NULL)
        if (is.null(fit) || is.null(fit$loglik) || isTRUE(fit$nevent == 0)) -Inf else fit$loglik[2]
      }
    },
    -Inf
  )

  if (!is.null(cache)) cache[[key]] <- stat_value
  return(stat_value)
}

#' Internal helper: Wrapper for the rgenoud Genetic Algorithm over Regularized Space
#'
#' @description
#' Sets up and runs `rgenoud::genoud` to optimize discrete indices over a regularized pool.
#'
#' @inheritParams find_cutpoint
#' @inheritParams find_cutpoint_number
#' @param target Continuous predictor vector.
#' @param numcut Number of cut-points.
#' @param time Survival time vector.
#' @param censor Survival event vector.
#' @param confound Optional covariate data frame.
#' @param gap Minimum distance between cut-points.
#' @param print.level Console output level.
#'
#' @return A list containing decoded `par` values and final fitness `value`.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G2.0} `requireNamespace("rgenoud")` check.
#' @srrstats {G5.8} Edge cases (zero gap, infinite range) return `NULL`.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom stats quantile as.formula na.omit
#' @importFrom survival coxph Surv
#' @noRd
.run_genetic_search <- function(target, numcut, time, censor, confound, nmin,
                                criterion, max.generations = 30, pop.size = 100,
                                gap = NULL, print.level = 0, candidate_cuts = NULL, ...) {
  if (!requireNamespace("rgenoud", quietly = TRUE)) {
    stop("The 'rgenoud' package is required. Please install it.", call. = FALSE)
  }

  # Establish the discrete grid coordinate matrix
  if (!is.null(candidate_cuts)) {
    grid_pool <- candidate_cuts
  } else {
    grid_probs <- seq(0.01, 0.99, by = 0.01)
    grid_pool <- sort(unique(stats::quantile(target, probs = grid_probs, na.rm = TRUE)))
  }

  G <- length(grid_pool)
  if (G < numcut) {
    return(NULL)
  }

  num_confound_vars <- if (is.null(confound)) 0 else ncol(confound)
  optimising_betas <- (criterion == "loglik")
  nvars <- if (optimising_betas) numcut * 2 + num_confound_vars else numcut

  # Map domains directly to valid integer boundaries (1 to G)
  domain_cuts <- matrix(rep(c(1, G), numcut), ncol = 2, byrow = TRUE)
  domain <- if (optimising_betas) rbind(domain_cuts, matrix(rep(c(-5, 5), nvars - numcut), ncol = 2, byrow = TRUE)) else domain_cuts

  # Compute starting index points evenly spaced across the matrix length
  initial_indices <- round(seq(2, G - 1, length.out = numcut))
  initial_values <- if (optimising_betas) c(initial_indices, rep(0, nvars - numcut)) else initial_indices

  if (is.null(gap)) {
    gap <- stats::quantile(sort(diff(sort(unique(stats::na.omit(target))))), probs = 0.5, na.rm = TRUE)
    if (is.na(gap) || gap == 0) gap <- 1e-4
  }

  loglik0 <- NA_real_
  if (criterion == "p_value") {
    null_formula_str <- "survival::Surv(time, censor) ~ 1"
    data_for_null_fit <- data.frame(time = time, censor = censor)
    if (num_confound_vars > 0) {
      colnames(confound) <- make.names(colnames(confound))
      data_for_null_fit <- cbind(data_for_null_fit, confound)
      null_formula_str <- paste(null_formula_str, "+", paste(colnames(confound), collapse = " + "))
    }
    null_fit <- tryCatch(survival::coxph(stats::as.formula(null_formula_str), data = data_for_null_fit), error = function(e) NULL)
    if (!is.null(null_fit)) loglik0 <- null_fit$loglik[2]
  }

  eval_cache <- new.env(hash = TRUE, parent = emptyenv())
  base_df <- data.frame(time = time, censor = censor)
  formula_str <- "survival::Surv(time, censor) ~ cut_design"
  if (num_confound_vars > 0) {
    colnames(confound) <- make.names(colnames(confound))
    base_df <- cbind(base_df, confound)
    formula_str <- paste(formula_str, "+", paste(colnames(confound), collapse = " + "))
  }

  optim_result <- tryCatch(rgenoud::genoud(
    fn = .obj, nvars = nvars, max = TRUE, pop.size = pop.size,
    max.generations = max.generations, wait.generations = 5, hard.generation.limit = TRUE,
    starting.values = initial_values, Domains = domain, print.level = print.level,
    data.type = 1, # DISCRETE INTEGER SEARCH MODE
    P9 = 0, # Turn off continuous local gradient optimization
    boundary.enforcement = 2,
    gradient.check = FALSE, # Bypass expensive fractional checking loops
    time = time, censor = censor, target = target, confound = confound,
    numcut = numcut, gap = gap, nmin = nmin, criterion = criterion,
    loglik0 = loglik0, cache = eval_cache, base_df = base_df,
    precompiled_formula = stats::as.formula(formula_str), grid_pool = grid_pool
  ), error = function(e) NULL)

  if (is.null(optim_result)) {
    return(NULL)
  }

  # Decode winning integer parameters back to real clinical coordinates
  final_indices <- round(optim_result$par[1:numcut])
  decoded_cuts <- sort(grid_pool[final_indices])

  return(list(
    par = decoded_cuts,
    value = optim_result$value
  ))
}

#' Internal helper: Genetic Model Selection for finding `max_cuts`
#'
#' @description
#' Evaluates models from 1 to `max_cuts` using the `rgenoud` evolutionary engine,
#' computing Information Criteria (AIC/BIC) to penalize complexity.
#'
#' @inheritParams find_cutpoint
#' @inheritParams find_cutpoint_number
#' @param userdata Cleaned survival data frame.
#'
#' @return A data.frame of model selection results.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.11} Helper for genetic-algorithm-based model selection.
#' @srrstats {G1.4a} Internal use only (`@noRd`).
#'
#' @importFrom survival coxph Surv
#' @importFrom stats as.formula
#' @importFrom cli cli_inform cli_alert_info
#' @noRd
.genetic_search_num <- function(userdata, max_cuts, nmin, criterion, covariates, max.generations, pop.size, ...) {
  n <- nrow(userdata)
  cov_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""
  base_formula_str <- paste("survival::Surv(time, event) ~ factor", cov_part)

  ic0 <- tryCatch(
    {
      fit0 <- survival::coxph(stats::as.formula(base_formula_str), data = userdata)
      .calc_ic(fit0, k = 1 + length(covariates), n = n, criterion = criterion)
    },
    error = function(e) {
      if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("Could not calculate IC for base model (0 cuts):", e$message))
      NA_real_
    }
  )

  results <- data.frame(num_cuts = 0, IC = ic0)
  results$cuts <- I(list(NULL))

  for (k_cuts in 1:max_cuts) {
    if (requireNamespace("cli", quietly = TRUE)) cli::cli_alert_info(paste("Running discrete genetic algorithm for", k_cuts, "cut-point(s)..."))

    ga_result <- tryCatch(
      {
        .run_genetic_search(
          target = userdata$factor, numcut = k_cuts, time = userdata$time, censor = userdata$event,
          confound = if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL,
          nmin = nmin, criterion = "loglik", max.generations = max.generations, pop.size = pop.size, ...
        )
      },
      error = function(e) {
        if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("Genetic algorithm failed for", k_cuts, "cut(s):", e$message))
        NULL
      }
    )

    if (!is.null(ga_result) && is.finite(ga_result$value) && ga_result$value > -.Machine$double.xmax) {
      ic_val <- .calc_ic(list(loglik = c(NA, ga_result$value)), k = k_cuts + length(covariates), n = n, criterion = criterion)
      results <- rbind(results, data.frame(num_cuts = k_cuts, IC = ic_val, cuts = I(list(sort(ga_result$par[1:k_cuts])))))
    } else {
      if (requireNamespace("cli", quietly = TRUE)) cli::cli_inform(paste("No valid cut-points found for", k_cuts, "cut(s) due to genetic algorithm failure or constraints."))
      results <- rbind(results, data.frame(num_cuts = k_cuts, IC = NA_real_, cuts = I(list(NULL))))
    }
  }
  names(results)[2] <- criterion
  return(results)
}

#' Internal helper: Calculate Information Criterion (AIC, AICc, or BIC)
#'
#' @description
#' Computes AIC, AICc, or BIC from a `coxph` object or log-likelihood.
#'
#' @inheritParams find_cutpoint_number
#' @param model Fitted `coxph` object or list with `loglik` component.
#' @param k Number of parameters (degrees of freedom).
#' @param n Sample size.
#'
#' @return Single numeric IC value (or `NA` on failure).
#'
#' @importFrom stats as.formula
#' @noRd
.calc_ic <- function(model, k, n, criterion) {
  if (is.null(model) || !is.list(model) || is.null(model$loglik)) {
    return(NA_real_)
  }
  logL <- model$loglik[2]
  if (is.na(logL)) {
    return(NA_real_)
  }
  if (criterion == "BIC" && (is.na(n) || n <= 0)) {
    return(NA_real_)
  }
  if (criterion == "BIC") {
    return(-2 * logL + k * log(n))
  } else if (criterion == "AICc") {
    if ((n - k - 1) <= 0) {
      return(NA_real_)
    }
    aic <- -2 * logL + 2 * k
    aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
    return(aicc)
  } else {
    return(-2 * logL + 2 * k)
  }
}
