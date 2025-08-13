#' Generalized Objective Function for Genetic Algorithm
#'
#' @description
#' This function is called by rgenoud::genoud. It calculates a fitness value
#' (e.g., log-likelihood, log-rank statistic) for a given set of cut-points.
#' rgenoud always maximizes, so all criteria are framed as maximization problems.
#'
#' @param params A numeric vector. The first `numcut` values are the cut-points.
#'   Subsequent values can be initial beta coefficients for the Cox model.
#' @param time,censor,target,confound The survival data and predictor variables.
#' @param numcut The number of cut-points.
#' @param gap The minimum required distance between cut-points.
#' @param nmin The minimum number of observations required in each group.
#' @param criterion The statistic to maximize: "loglik", "logrank", "p_value", or "hazard_ratio".
#'
#' @return A single numeric value representing the fitness of the solution.
#' @noRd
.obj <- function(params, time, censor, target, confound, numcut, gap, nmin, criterion) {
  # --- 1. Extract cut-points and perform validation ---
  cutoff <- params[1:numcut]

  # Safeguards to ensure valid grouping
  if (numcut > 1 && min(diff(sort(cutoff))) < gap) return(-Inf) # Cuts are too close
  cut_design <- factor(cut(target, breaks = c(-Inf, sort(cutoff), Inf), labels = FALSE, include.lowest = TRUE))
  if (length(table(cut_design)) != (numcut + 1)) return(-Inf) # Incorrect number of groups formed
  if (min(table(cut_design)) < nmin) return(-Inf) # A group is too small

  # --- 2. Prepare data and formula for model fitting ---
  data_for_fit <- data.frame(time = time, censor = censor, cut_design = cut_design)
  if (!is.null(confound) && ncol(confound) > 0) {
    # Sanitize covariate names to be valid R variable names
    colnames(confound) <- make.names(colnames(confound))
    data_for_fit <- cbind(data_for_fit, confound)
  }

  formula_str <- "survival::Surv(time, censor) ~ cut_design"
  if (!is.null(confound) && ncol(confound) > 0) {
    formula_str <- paste(formula_str, "+", paste(colnames(confound), collapse = " + "))
  }
  fit_formula <- as.formula(formula_str)

  # --- 3. Calculate the fitness value based on the chosen criterion ---
  stat_value <- switch(
    criterion,
    "logrank" = {
      fit <- tryCatch(survival::survdiff(fit_formula, data = data_for_fit), error = function(e) NULL)
      if (is.null(fit)) -Inf else fit$chisq
    },
    "hazard_ratio" = {
      # This is only valid for numcut = 1. A check in the main find_cutpoint function prevents misuse.
      fit <- tryCatch(survival::coxph(fit_formula, data = data_for_fit), error = function(e) NULL)
      if (is.null(fit)) -Inf else summary(fit)$conf.int[1, "exp(coef)"]
    },
    "p_value" = {
      fit <- tryCatch(survival::coxph(fit_formula, data = data_for_fit), error = function(e) NULL)
      # We maximize (1 - p_value), which is equivalent to minimizing p_value.
      if (is.null(fit) || is.null(summary(fit)$logtest)) -Inf else 1 - summary(fit)$logtest["pvalue"]
    },
    "loglik" = {
      # This efficient path is for find_cutpoint_number. It uses initial betas to avoid a full fit.
      beta_init <- params[(numcut + 1):length(params)]
      fit <- tryCatch(survival::coxph(fit_formula, data = data_for_fit, init = beta_init, iter.max = 0), error = function(e) NULL)
      if (is.null(fit) || is.null(fit$loglik)) -Inf else fit$loglik[2]
    },
    # Default case if criterion is unknown, returns a very poor fitness score
    -Inf
  )

  return(stat_value)
}


#' Wrapper for the rgenoud Genetic Algorithm
#'
#' @description
#' Sets up and runs the genetic algorithm using rgenoud::genoud to find the
#' optimal cut-points for a given criterion.
#'
#' @param target The continuous predictor variable vector.
#' @param numcut The number of cut-points to find.
#' @param criterion The optimization criterion.
#' @param ... Other parameters passed to .obj and genoud.
#'
#' @return A `genoud` object containing the results of the optimization.
#' @noRd
.run_genetic_search <- function(target, numcut, time, censor, confound, nmin, criterion,
                                numgen = 15, gap = NULL, print.level = 0, ...) {
  # --- 1. Determine number of variables to optimize ---
  num_confound_vars <- if (is.null(confound) || all(is.na(confound))) 0 else ncol(confound)

  # For 'loglik', we optimize both cuts and model betas. For others, just the cuts.
  optimizing_betas <- (criterion == "loglik")
  nvars <- if (optimizing_betas) {
    numcut + numcut + num_confound_vars # cuts + group betas + confounder betas
  } else {
    numcut # Just the cut-points
  }

  # --- 2. Define search space (domains) ---
  domain_cuts <- matrix(rep(range(target, na.rm=TRUE), numcut), ncol = 2, byrow = TRUE)
  if(any(is.infinite(domain_cuts))) return(NULL)

  domain <- if (optimizing_betas) {
    # Define a wide search space for beta coefficients
    domain_betas <- matrix(rep(c(-5, 5), nvars - numcut), ncol = 2, byrow = TRUE)
    rbind(domain_cuts, domain_betas)
  } else {
    domain_cuts
  }

  # --- 3. Set sensible starting values ---
  initial_cuts <- stats::quantile(target, probs = seq(0, 1, length.out = numcut + 2), na.rm = TRUE)[2:(numcut + 1)]
  initial_values <- if (optimizing_betas) {
    initial_betas <- rep(0, nvars - numcut)
    c(initial_cuts, initial_betas)
  } else {
    initial_cuts
  }
  if(any(is.na(initial_values))) return(NULL)

  # --- 4. Define a minimum gap between cut-points ---
  if (is.null(gap)) {
    # Use the median difference between sorted unique values as a default gap
    gap <- stats::quantile(sort(diff(sort(unique(na.omit(target))))), probs = 0.5, na.rm = TRUE)
    if(is.na(gap) || gap == 0) gap <- 1e-4 # Fallback for edge cases
  }

  # --- 5. Run the Genetic Algorithm ---
  if (!requireNamespace("rgenoud", quietly = TRUE)) {
    stop("The 'rgenoud' package is required for this function. Please install it.", call. = FALSE)
  }

  optim_result <- tryCatch(rgenoud::genoud(
    fn = .obj, nvars = nvars, max = TRUE, pop.size = 100,
    max.generations = numgen,
    wait.generations = 10, hard.generation.limit = TRUE,
    starting.values = initial_values, Domains = domain, print.level = print.level,
    # Pass all necessary arguments to the objective function .obj
    time = time, censor = censor, target = target, confound = confound,
    numcut = numcut, gap = gap, nmin = nmin, criterion = criterion
  ), error = function(e) NULL)

  return(optim_result)
}

#' Calculate Information Criterion (AIC, AICc, or BIC) for Survival Models
#'
#' @description
#' Calculates AIC, BIC, or AICc from a fitted model object or log-likelihood value.
#'
#' @param model A fitted `coxph` model object, or a list containing the log-likelihood.
#' @param k The number of parameters in the model (degrees of freedom).
#' @param n The sample size.
#' @param criterion The information criterion to calculate.
#'
#' @return A single numeric value for the requested information criterion.
#' @noRd
.calc_ic <- function(model, k, n, criterion) {
  # This function can now accept a raw logL value via a list,
  # which is useful when we get it directly from the GA optimizer.
  if(is.null(model) || !is.list(model)) return(NA_real_)
  logL <- model$loglik[2]
  if(is.null(logL) || is.na(logL)) return(NA_real_)

  if (criterion == "BIC") {
    return(-2 * logL + k * log(n))
  } else if (criterion == "AICc") {
    # AICc requires n > k + 1
    if ((n - k - 1) <= 0) return(NA_real_)
    aic <- -2 * logL + 2 * k
    aicc <- aic + (2 * k * (k + 1)) / (n - k - 1)
    return(aicc)
  } else { # AIC
    return(-2 * logL + 2 * k)
  }
}

#' Infix operator for providing a default value for NULL.
#' @param a The value to check.
#' @param b The default value to use if a is NULL.
#' @name or-operator
#' @export
`%||%` <- function(a, b) if (is.null(a)) b else a
