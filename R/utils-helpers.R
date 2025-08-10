#' Objective Function for Genetic Algorithm
#'
#' @noRd
.obj <- function(params, time, censor, target, confound, numcut, gap, nmin) {
  cutoff <- params[1:numcut]

  # Safeguards
  if (numcut > 1 && min(diff(sort(cutoff))) < gap) return(-Inf)
  cut_design <- factor(cut(target, breaks = c(-Inf, sort(cutoff), Inf), labels = FALSE, include.lowest = TRUE))
  if (length(table(cut_design)) != (numcut + 1)) return(-Inf)
  if (min(table(cut_design)) < nmin) return(-Inf)

  beta_init <- params[(numcut + 1):length(params)]

  if (is.null(confound) || ncol(confound) == 0) {
    fit <- tryCatch(survival::coxph(survival::Surv(time, censor) ~ cut_design, init = beta_init, iter.max = 0), error = function(e) NULL)
  } else {
    fit <- tryCatch(survival::coxph(survival::Surv(time, censor) ~ cut_design + confound, init = beta_init, iter.max = 0), error = function(e) NULL)
  }

  if(is.null(fit)) return(-Inf)
  return(fit$loglik[2])
}


#' Wrapper for the rgenoud Genetic Algorithm
#'
#' @noRd
.maxloglik <- function(target, numcut, time, censor, confound, nmin, numgen = 15, gap = NULL, print.level = 0) {
  if(length(na.omit(target)) < 20) return(NULL)

  num_confound_vars <- if (is.null(confound) || all(is.na(confound))) 0 else ncol(confound)
  nvars <- numcut + numcut + num_confound_vars

  domain_cuts <- matrix(rep(range(target, na.rm=TRUE), numcut), ncol = 2, byrow = TRUE)
  if(any(is.infinite(domain_cuts))) return(NULL)

  domain_betas <- matrix(rep(c(-5, 5), nvars - numcut), ncol = 2, byrow = TRUE)
  domain <- rbind(domain_cuts, domain_betas)

  initial_cuts <- stats::quantile(target, probs = seq(0, 1, length.out = numcut + 2), na.rm = TRUE)[2:(numcut + 1)]
  initial_betas <- rep(0, nvars - numcut)
  initial_values <- c(initial_cuts, initial_betas)

  if(any(is.na(initial_values))) return(NULL)

  if (is.null(gap)) {
    gap <- stats::quantile(sort(diff(sort(na.omit(target)))), probs = 0.5, na.rm = TRUE)
    if(is.na(gap) || gap == 0) gap <- 1e-4
  }

  optim_result <- tryCatch(rgenoud::genoud(
    fn = .obj, nvars = nvars, max = TRUE, pop.size = 100,
    max.generations = numgen,
    wait.generations = 10, hard.generation.limit = TRUE,
    starting.values = initial_values, Domains = domain, print.level = print.level,
    time = time, censor = censor, target = target, confound = confound,
    numcut = numcut, gap = gap, nmin = nmin
  ), error = function(e) NULL)

  return(optim_result)
}


#' Calculate Information Criterion (AIC or BIC)
#'
#' @noRd
.calc_ic <- function(fit, k, n_obs, type) {
    logLik_val <- if(inherits(fit, "coxph")) fit$loglik[2] else as.numeric(stats::logLik(fit))
    if (type == "AIC") return(-2 * logLik_val + 2 * k)
    if (type == "BIC") return(-2 * logLik_val + log(n_obs) * k)
}