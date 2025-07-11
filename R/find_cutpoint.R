#' Find an Optimal Cut-point for a Continuous Predictor
#'
#' @description
#' Determines the optimal cut-point(s) for a continuous variable. It can use
#' a "systematic" search or a "genetic" algorithm (survival only). The analysis
#' can be adjusted for specified covariates.
#'
#' @param data A data frame containing the variables.
#' @param predictor The name of the predictor variable, as a string.
#' @param covariates A character vector of covariate names to adjust for in the model.
#' @param outcome_time The name of the time variable, as a string (for survival).
#' @param outcome_event The name of the event variable, as a string (for survival).
#' @param outcome_binary The name of the binary outcome variable, as a string.
#' @param method The algorithm to use: "systematic" or "genetic".
#' @param num_cuts The number of cut-points to find.
#' @param nmin The minimum number of observations in each group.
#' @param vote_methods For the systematic method, a character vector specifying
#'   which criteria to include in the majority vote.
#' @param use_parallel Logical. If TRUE, uses multiple CPU cores.
#' @param ... Additional arguments passed to the genetic algorithm.
#'
#' @return An object containing the results of the cut-point analysis.
#' @importFrom foreach %dopar%
#' @importFrom stats na.omit as.formula pchisq
#' @importFrom graphics points
#' @export
find_cutpoint <- function(data, predictor, covariates = NULL,
                          outcome_time = NULL, outcome_event = NULL,
                          outcome_binary = NULL,
                          method = "systematic",
                          num_cuts = 1, nmin = 20,
                          vote_methods = NULL,
                          use_parallel = FALSE, ...) {

  # --- 1. Input Validation and Data Prep ---
  method <- match.arg(method, choices = c("systematic", "genetic"))
  if (is.null(predictor)) stop("A 'predictor' variable must be specified as a string.", call. = FALSE)

  if (!is.null(outcome_time) && !is.null(outcome_event)) {
    analysis_type <- "survival"
    required_vars <- c(predictor, outcome_time, outcome_event, covariates)
    if(is.null(vote_methods)) vote_methods <- c("logrank", "cox", "hr")
  } else if (!is.null(outcome_binary)) {
    analysis_type <- "logistic"
    required_vars <- c(predictor, outcome_binary, covariates)
    if(is.null(vote_methods)) vote_methods <- c("youden", "auc", "closest_topleft")
  } else {
    stop("You must specify outcome variables for this method.", call. = FALSE)
  }
  if (!all(required_vars %in% names(data))) {
    missing_cols <- required_vars[!required_vars %in% names(data)]
    stop(paste("The following specified columns were not found:", paste(missing_cols, collapse = ", ")), call. = FALSE)
  }

  all_vars <- c(predictor, covariates, outcome_time, outcome_event, outcome_binary)
  userdata <- data[, unique(all_vars), drop = FALSE]
  userdata <- na.omit(userdata)

  names(userdata)[names(userdata) == predictor] <- "factor"
  if(analysis_type == "survival") {
    names(userdata)[names(userdata) == outcome_time] <- "time"
    names(userdata)[names(userdata) == outcome_event] <- "event"
  } else {
    names(userdata)[names(userdata) == outcome_binary] <- "outcome"
  }

  cli::cli_alert_info("Starting cut-point analysis: method = {.strong {method}}")
  if (!is.null(covariates)) cli::cli_alert_info("Adjusting for covariates: {paste(covariates, collapse=', ')}")

  n_total <- nrow(userdata)
  if (nmin < 1 && nmin > 0) {
    nmin_abs <- ceiling(nmin * n_total)
    cli::cli_alert_info("Interpreting nmin = {nmin} as a proportion. Minimum group size set to {nmin_abs}.")
    nmin <- nmin_abs
  } else if (nmin >= 1) {
    nmin <- as.integer(nmin)
  } else {
    stop("'nmin' must be a positive number.", call. = FALSE)
  }

  # --- 2. Route to the correct method ---
  if (method == "systematic") {

    if (num_cuts != 1) stop("The 'systematic' method only supports num_cuts = 1.", call. = FALSE)
    userdata <- userdata[order(userdata$factor), ]
    n <- nrow(userdata)

    covariate_formula_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""

    .eval_survival_cut <- function(cut_val, data_in, nmin_val) {
      data_in$group <- factor(ifelse(data_in$factor <= cut_val, 1, 2))
      if (length(unique(data_in$group)) < 2 || min(table(data_in$group)) < nmin_val) return(NULL)

      formula_str <- paste("survival::Surv(time, event) ~ group", covariate_formula_part)
      fit <- tryCatch(survival::coxph(as.formula(formula_str), data = data_in),
                      warning = function(w) NULL, error = function(e) NULL)

      if (is.null(fit)) return(NULL)

      model_summary <- summary(fit)
      if (is.null(model_summary$conf.int)) return(NULL)

      lr_test <- survival::survdiff(survival::Surv(time, event) ~ group, data = data_in)

      data.frame(Cut1 = cut_val, LogRank_p = 1 - pchisq(lr_test$chisq, df = 1),
                 HR = model_summary$conf.int["group2", "exp(coef)"],
                 Cox_p = summary(fit)$sctest["pvalue"])
    }
    .eval_logistic_cut <- function(cut_val, data_in, nmin_val) {
      data_in$group <- factor(ifelse(data_in$factor <= cut_val, 1, 2))
      if (length(unique(data_in$group)) < 2 || min(table(data_in$group)) < nmin_val) return(NULL)

      formula_str <- paste("outcome ~ group", covariate_formula_part)
      fit <- stats::glm(as.formula(formula_str), data = data_in, family = stats::binomial())

      predicted_probs <- stats::predict(fit, type = "response")
      roc_obj <- suppressMessages(pROC::roc(data_in$outcome, predicted_probs, quiet = TRUE))

      coords <- pROC::coords(roc_obj, x = "best", best.method = "youden", ret = c("threshold", "sensitivity", "specificity"))

      data.frame(Cut1 = cut_val,
                 AUC = as.numeric(pROC::auc(roc_obj)),
                 Youden = coords$sensitivity + coords$specificity - 1,
                 Closest_Topleft = sqrt((1 - coords$sensitivity)^2 + (coords$specificity)^2),
                 Sens_Spec_Product = coords$sensitivity * coords$specificity)
    }

    search_grid <- unique(userdata$factor[nmin:(n - nmin)])
    eval_func <- if (analysis_type == "survival") .eval_survival_cut else .eval_logistic_cut

    if (use_parallel) {
      cores <- parallel::detectCores()
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      cli::cli_alert_info("Using {cores} cores for parallel processing...")
    } else {
      foreach::registerDoSEQ()
    }

    allcut <- foreach::foreach(cut_val = search_grid, .combine = 'rbind') %dopar% {
      eval_func(cut_val, data_in = userdata, nmin_val = nmin)
    }

    if (is.null(allcut) || nrow(allcut) == 0) stop("No valid cut-points found.", call. = FALSE)

    if (analysis_type == "survival") {
      all_votes <- list()
      if("logrank" %in% vote_methods) all_votes$logrank <- allcut$Cut1[which.min(allcut$LogRank_p)]
      if("cox" %in% vote_methods) all_votes$cox <- allcut$Cut1[which.min(allcut$Cox_p)]
      if("hr" %in% vote_methods) all_votes$hr <- allcut$Cut1[which.max(abs(log(allcut$HR)))]

      if("maxstat" %in% vote_methods) {
        if (!requireNamespace("maxstat", quietly = TRUE)) stop("Package 'maxstat' is required for this method.")
        ms_fit <- maxstat::maxstat.test(survival::Surv(time, event) ~ factor, data=userdata, smethod="LogRank", pmethod="none")
        all_votes$maxstat <- ms_fit$estimate
      }

      best_cuts_summary <- data.frame(Criterion = names(all_votes), Cut_Point = unlist(all_votes))

    } else {
      all_votes <- list()
      if("youden" %in% vote_methods) all_votes$youden <- allcut$Cut1[which.max(allcut$Youden)]
      if("auc" %in% vote_methods) all_votes$auc <- allcut$Cut1[which.max(allcut$AUC)]
      if("closest_topleft" %in% vote_methods) all_votes$closest_topleft <- allcut$Cut1[which.min(allcut$Closest_Topleft)]
      if("sens_spec_product" %in% vote_methods) all_votes$sens_spec_product <- allcut$Cut1[which.max(allcut$Sens_Spec_Product)]

      if ("distribution" %in% vote_methods) {
        dist_data <- data.frame(factor = data[[predictor]])
        dist_res <- tryCatch(.distribution_cut(dist_data), error = function(e) NULL)
        if (!is.null(dist_res)) all_votes$distribution <- dist_res$optimal_cut
      }

      selected_method_cuts <- unlist(all_votes[names(all_votes) %in% vote_methods])
      best_cuts_summary <- data.frame(Criterion = names(selected_method_cuts), Cut_Point = selected_method_cuts)
    }

    vote_table <- table(best_cuts_summary$Cut_Point)
    best_by_vote <- as.numeric(names(vote_table)[which.max(vote_table)])

    output <- list(allcut = allcut, best_cuts_summary = best_cuts_summary, best_by_vote = best_by_vote,
                   parameters = list(analysis_type = analysis_type, nmin = nmin, num_cuts = num_cuts, covariates = covariates),
                   userdata = userdata)
    class(output) <- "find_cutpoint_systematic"
    cli::cli_alert_success("Systematic search complete.")
    return(output)

  } else if (method == "genetic") {

    if (analysis_type != "survival") stop("The 'genetic' method is only available for survival outcomes.", call. = FALSE)
    if (!exists(".maxloglik", mode = "function")) stop("Helper function .maxloglik not found.", call. = FALSE)

    confound_data <- if (!is.null(covariates)) as.matrix(userdata[, covariates]) else NULL

    gen_result <- .maxloglik(
      target = userdata$factor, numcut = num_cuts, time = userdata$time,
      censor = userdata$event, confound = confound_data, nmin = nmin, ...
    )
    if (is.null(gen_result)) stop("Optimization failed, possibly due to sparse data.", call. = FALSE)

    final_cuts <- sort(gen_result$par[1:num_cuts])
    userdata$group <- factor(cut(userdata$factor, breaks = c(-Inf, final_cuts, Inf)), labels = paste0("G", 1:(num_cuts+1)))

    covariate_formula_part <- if (!is.null(covariates)) paste(" +", paste(covariates, collapse = " + ")) else ""
    final_formula <- as.formula(paste("survival::Surv(time, event) ~ group", covariate_formula_part))
    fit <- survival::coxph(final_formula, data = userdata)

    output <- list(summary = summary(fit), optimal_cuts = final_cuts, optim_details = gen_result,
                   userdata = userdata, parameters = list(analysis_type = analysis_type, num_cuts = num_cuts, covariates = covariates))
    class(output) <- "find_cutpoint_genetic"
    cli::cli_alert_success("Genetic algorithm complete.")
    return(output)
  }
}

# --- Internal Helper Function for Distribution Method ---
.distribution_cut <- function(userdata) {
  if (!requireNamespace("flexmix", quietly = TRUE)) {
    stop("Package 'flexmix' is required for the distribution method.", call. = FALSE)
  }
  if (length(unique(userdata$factor)) < 3) return(NULL)

  fit <- tryCatch(flexmix::flexmix(factor ~ 1, data = userdata, k = 2), error = function(e) NULL)
  if (is.null(fit)) return(NULL)

  params <- flexmix::parameters(fit)
  priors <- flexmix::prior(fit)
  f <- function(x) {
    (priors[1] * stats::dnorm(x, mean = params[1, 1], sd = params[2, 1])) -
      (priors[2] * stats::dnorm(x, mean = params[1, 2], sd = params[2, 2]))
  }
  intersect <- tryCatch(stats::uniroot(f, lower = min(params[1,]), upper = max(params[1,]))$root,
                        error = function(e) NA)

  return(list(optimal_cut = intersect, model = fit))
}
