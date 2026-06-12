# ===================================================================
# MAIN FUNCTION: FIND OPTIMAL CUT-POINTS FOR SURVIVAL DATA
# ===================================================================

#' Find Optimal Cut-points for Survival Data
#'
#' @description
#' Finds optimal cut-point(s) for a continuous predictor in a
#' time-to-event (survival) analysis. Uses systematic search (1–2
#' cuts) or a genetic algorithm (any number of cuts). Features high-speed
#' integer partitioning via compiled C++ vector assignments and automated
#' quantile grid downsampling.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G1.1} Implements genetic + systematic search for optimal
#'   multi-cutpoint survival groupings.
#' @srrstats {G1.0} References provided for Cox, log-rank, genetic optimisation.
#' @srrstats {G1.3} Systematic grid search (1–2 cuts) and `rgenoud` global
#'   optimisation documented.
#' @srrstats {G1.5} Compared with `cutpointr` and `survminer` in package
#'   vignette.
#' @srrstats {G1.6} Numerical stability via `survival::coxph` and `rgenoud`;
#'   edge cases return `NA`.
#' @srrstats {G2.3a} Uses `match.arg()` to validate `method` and `criterion`
#'   arguments.
#' @srrstats {G2.3b} Uses `as.formula` and safe subsetting (NSE safe).
#' @srrstats {G2.4} `NA` removed via `stats::na.omit()`.
#' @srrstats {G2.5} Factor ordering is handled by `cut()` which creates ordered
#'   factors by default.
#' @srrstats {RE4.0} Returns class `find_cutpoint` with model details.
#' @srrstats {G5.2} `optimal_cuts` and `optimal_stat` are `NA` when no valid
#'   solution found (Graceful failure).
#' @srrstats {G2.6} Validates inputs via helpers.
#' @srrstats {G2.8} Informative errors via `cli::cli_abort()`.
#' @srrstats {G5.2} Warnings via `cli::cli_alert_warning()`.
#' @srrstats {G5.2} Graceful degradation via `na_result()` for empty data or
#'   model failures.
#' @srrstats {G2.13} `cli_abort()` for invalid input (missingness checks).
#' @srrstats {G2.15} Explicit checks prevent passing missing data to analytic
#'   functions (via `na.omit`).
#' @srrstats {G2.0} `cli_abort()` checks for missing `rgenoud` dependency.
#' @srrstats {G2.14c} `NA` propagation controlled.
#' @srrstats {RE6.0} `plot()` method provided.
#' @srrstats {RE6.1} `plot()` method is an S3 generic dispatch.
#' @srrstats {G1.4} All parameters and return values documented.
#' @srrstats {RE4.2} Model selection via log-rank, HR, or p-value.
#' @srrstats {RE6.2} Visualises fitted values (survival curves).
#' @srrstats {RE1.0} Implements optimal cut-point algorithm.
#' @srrstats {RE1.1} Assumes PH; check `summary()` for `cox.zph`.
#' @srrstats {RE1.4} Cox PH assumption test via `summary(fit)$cox_zph`.
#' @srrstats {RE2.0} Transformations documented in details.
#' @srrstats {RE2.1} Estimates/SEs from `coxph` in `summary()`.
#' @srrstats {RE3.0} `tryCatch` checks model convergence (convergence warnings).
#' @srrstats {RE2.4} Collinearity checks via model constraints.
#' @srrstats {RE2.4a} Checks for collinearity among predictors.
#' @srrstats {RE2.4b} Checks for collinearity between X and Y.
#'
#' @details
#' `method = "systematic"`: grid search respecting `nmin`. Optimised via internal quantiles.
#' `method = "genetic"`: `rgenoud` global optimisation.
#' Systematic search is slow for `num_cuts > 2`; use `genetic`.
#' Core vector partitions are calculated in compiled C++ via `Rcpp` for optimal performance.
#'
#' @references
#' Altman, D. G., Lausen, B., Sauerbrei, W., & Schumacher,
#' M. (1994). Dangers of Using “Optimal” Cutpoints in the Evaluation of
#' Prognostic Factors. *JNCI: Journal of the National Cancer Institute*,
#' 86(11), 829–835. \doi{10.1093/jnci/86.11.829}
#'
#' Cox, D. R. (1972). Regression Models and Life-Tables. *Journal
#' of the Royal Statistical Society: Series B (Methodological)*, 34(2),
#' 187–202. \doi{10.1111/j.2517-6161.1972.tb00899.x}
#'
#' Mantel, N. (1966). Evaluation of survival data and two new
#' rank order statistics arising in its consideration. *Cancer
#' Chemotherapy Reports*, 50(3).
#'
#' Mebane Jr, W. R., & Sekhon, J. S. (2011). Genetic
#' Optimization Using Derivatives: The rgenoud Package for R.
#' *Journal of Statistical Software*, 42, 1–26.
#' \doi{10.18637/jss.v042.i11}
#'
#' @param data A data frame containing the analysis variables.
#' @param predictor The continuous predictor variable name (character).
#' @param outcome_time The time-to-event variable name (character).
#' @param outcome_event The event status variable name (character, 0 or 1).
#' @param num_cuts The number of cut-points to find. Default is 1.
#' @param method Algorithm search type: `"systematic"` or `"genetic"`.
#' @param criterion The statistic to optimise: `"logrank"`, `"hazard_ratio"`, or `"p_value"`.
#' @param covariates Character vector of covariate names (optional).
#' @param nmin Min. group size (integer count or proportion).
#' @param seed Optional integer seed for reproducible genetic search.
#' @param max.generations Max generations for genetic algorithm. If `NULL`, dynamically scales.
#' @param pop.size Population size for genetic algorithm. If `NULL`, dynamically scales.
#' @param n_perm Number of permutations to run for an adjusted p-value. Default is 0.
#' @param n_cores Number of CPU cores for parallel permutations. Default is 1.
#' @param use_cpp Logical. Checks and calls compiled C++ routines via `Rcpp`. Default is `TRUE`.
#' @param grid_by Percentile step increment for systematic grid downsampling. Default is 0.01.
#' @param quiet Logical. If `TRUE`, suppresses operational console alerts.
#' @param candidate_cuts Optional vector of pre-filtered cuts defining a narrow search space.
#' @param ... Additional arguments passed directly to `rgenoud::genoud`.
#'
#' @return An object of class `find_cutpoint` containing the
#'      optimal cut-points, statistic, and analysis parameters.
#' @useDynLib OptSurvCutR, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @export
#'
#' @examples
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'
#'   # Create a lightweight, reproducible simulation baseline dataset
#'   set.seed(42)
#'   sim_data <- data.frame(
#'     time = rexp(30, rate = 0.1),
#'     event = sample(c(0, 1), 30, replace = TRUE),
#'     biomarker = rnorm(30, mean = 5, sd = 1.5)
#'   )
#'
#'   # Execute an exhaustive systematic threshold discovery sweep
#'   fit <- find_cutpoint(
#'     data = sim_data,
#'     predictor = "biomarker",
#'     outcome_time = "time",
#'     outcome_event = "event",
#'     num_cuts = 1,
#'     method = "systematic",
#'     criterion = "logrank",
#'     nmin = 5,
#'     quiet = TRUE
#'   )
#'   print(fit)
#' }
find_cutpoint <- function(data, predictor, outcome_time, outcome_event,
                          num_cuts = 1, method = c("systematic", "genetic"),
                          criterion = c("logrank", "hazard_ratio", "p_value"),
                          covariates = NULL, nmin = 20, seed = NULL,
                          max.generations = NULL, pop.size = NULL,
                          n_perm = 0, n_cores = 1, use_cpp = TRUE,
                          grid_by = 0.01, quiet = FALSE, candidate_cuts = NULL, ...) {
  method <- match.arg(method)
  criterion <- match.arg(criterion)

  .validate_find_cutpoint_inputs(data, predictor, outcome_time, outcome_event, num_cuts, method, criterion, covariates)
  userdata <- .prepare_cutpoint_data(data, predictor, outcome_time, outcome_event, covariates)

  # Graceful NA return helper
  na_return <- function(nmin_abs = nmin) {
    res <- list(
      optimal_cuts = rep(NA_real_, num_cuts), optimal_stat = NA_real_, all_stats = NULL, userdata = userdata,
      parameters = list(
        method = method, analysis_type = "survival", predictor = predictor, num_cuts = num_cuts,
        criterion = criterion, covariates = covariates, nmin = nmin_abs, max.generations = max.generations,
        pop.size = pop.size, use_cpp = use_cpp, grid_by = grid_by, quiet = quiet
      )
    )
    class(res) <- "find_cutpoint"
    return(res)
  }

  if (nrow(userdata) == 0) {
    if (!quiet) cli::cli_inform("No complete cases found after removing NAs.")
    return(na_return())
  }

  val_res <- .validate_data_conditions(userdata, nmin, num_cuts, outcome_event, quiet)
  if (!val_res$valid) {
    return(na_return(val_res$nmin_abs))
  }
  nmin_abs <- val_res$nmin_abs

  if (is.factor(userdata$factor) || is.character(userdata$factor)) {
    cli::cli_abort(c(
      "x" = "Unordered factors or characters are not allowed as continuous optimization inputs.",
      "i" = "Variable '{predictor}' must be passed as a continuous numeric vector."
    ))
  }

  if (use_cpp && !exists("cpp_get_group_assignments", mode = "function")) {
    if (!quiet) cli::cli_alert_warning("Compiled C++ binary not loaded. Falling back gracefully to native R vector processing.")
    use_cpp <- FALSE
  }

  if (!is.null(seed)) set.seed(seed)

  # Compile underlying search lattice to verify data node supply bounds
  search_grid_probs <- if (!is.null(grid_by)) seq(grid_by, 1 - grid_by, by = grid_by) else seq(0.01, 0.99, by = 0.01)
  base_grid <- if (!is.null(candidate_cuts)) candidate_cuts else sort(unique(stats::quantile(userdata$factor, probs = search_grid_probs, na.rm = TRUE)))

  # ===================================================================
  # UX COHORT HEADROOM CHECK (Dynamic Lot Defense)
  # ===================================================================
  n_total <- nrow(userdata)
  unique_values_needed <- (num_cuts + 1) * nmin_abs
  actual_available_nodes <- length(base_grid)

  if (n_total < unique_values_needed || actual_available_nodes < num_cuts) {
    cli::cli_abort(c(
      "x" = "Insufficient data density for {num_cuts} cut-point(s) at nmin = {nmin_abs}.",
      "i" = "Model requires at least {num_cuts} unique partition boundaries, but only {actual_available_nodes} survived duplicate filtering.",
      "*" = "Action: Reduce 'num_cuts' or lower your 'nmin' threshold."
    ))
  }

  # --- UX AUTO-SCALING FOR REGULARIZED SPACE ENGINE OVERRIDES ---
  if (is.null(pop.size) || pop.size == 100) {
    pop.size <- switch(as.character(num_cuts),
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
    max.generations <- switch(as.character(num_cuts),
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

  extra_args <- list(...)

  # Force hard boundary enforcement inside global genetic evaluations
  if (method == "genetic" && !("boundary.enforcement" %in% names(extra_args))) {
    extra_args$boundary.enforcement <- 2
  }

  if (method == "systematic") {
    real_res <- do.call(.systematic_search, c(
      list(
        userdata = userdata, num_cuts = num_cuts, criterion = criterion,
        covariates = covariates, nmin = nmin_abs, predictor_name = predictor,
        use_cpp = use_cpp, grid_by = grid_by, quiet = quiet, candidate_cuts = candidate_cuts
      ),
      extra_args
    ))
  } else {
    if (!quiet) cli::cli_alert_info("Starting regularized genetic search for {num_cuts} cut(s)...")
    confound_df <- if (!is.null(covariates)) userdata[, covariates, drop = FALSE] else NULL

    real_res <- do.call(.run_genetic_search, c(
      list(
        target = userdata$factor, numcut = num_cuts, time = userdata$time,
        censor = userdata$event, confound = confound_df, nmin = nmin_abs,
        criterion = criterion, max.generations = max.generations, pop.size = pop.size,
        use_cpp = use_cpp, print.level = 0, candidate_cuts = candidate_cuts
      ),
      extra_args
    ))

    if (is.null(real_res) || !is.finite(real_res$value) || real_res$value <= -.Machine$double.xmax) {
      if (!quiet) cli::cli_inform("Genetic algorithm found no valid solution.")
      real_res <- list(optimal_cuts = rep(NA_real_, num_cuts), optimal_stat = NA_real_, all_stats = NULL)
    } else {
      real_res <- list(optimal_cuts = sort(real_res$par[1:num_cuts]), optimal_stat = real_res$value, all_stats = NULL)
    }
  }

  p_perm <- NA
  if (n_perm > 0 && !anyNA(real_res$optimal_cuts)) {
    if (requireNamespace("cli", quietly = TRUE) && !quiet) cli::cli_alert_info("Running {n_perm} permutations to calculate adjusted p-value...")

    p_perm <- .run_permutations(
      time_vec = userdata$time,
      censor_vec = userdata$event,
      predictor_vec = userdata$factor,
      num_cuts = num_cuts,
      method = method,
      criterion = criterion,
      covariates = covariates,
      userdata_full = userdata,
      nmin = nmin_abs,
      obs_stat = real_res$optimal_stat,
      n_perm = n_perm,
      n_cores = n_cores,
      max.generations = max.generations,
      pop.size = pop.size,
      use_cpp = use_cpp,
      grid_by = grid_by,
      candidate_cuts = candidate_cuts
    )
  }

  output <- list(
    optimal_cuts = real_res$optimal_cuts,
    optimal_stat = real_res$optimal_stat,
    permuted_p_value = p_perm,
    n_perm = n_perm,
    all_stats = real_res$all_stats,
    userdata = userdata,
    parameters = list(
      method = method, analysis_type = "survival", predictor = predictor,
      num_cuts = num_cuts, criterion = criterion, covariates = covariates,
      nmin = nmin_abs, max.generations = max.generations, pop.size = pop.size,
      use_cpp = use_cpp, grid_by = grid_by, quiet = quiet
    )
  )

  class(output) <- "find_cutpoint"
  return(output)
}
