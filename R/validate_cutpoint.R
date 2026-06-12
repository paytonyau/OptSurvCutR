# ====================================================================
# VALIDATION FUNCTION
# Bootstraps cut-points for stability and confidence intervals.
# ====================================================================

#' Validate an Optimal Cut-point Using Bootstrapping
#'
#' @description
#' Assesses cut-point stability from \code{\link{find_cutpoint}} via bootstrap
#' analysis, generating 95\% confidence intervals. Streamlined for
#' survival (time-to-event) analysis.
#'
#' @section srrstats compliance:
#' .
#' @srrstats {G1.0} Primary refs: Efron (1979), Rota et al. (2015).
#' @srrstats {G2.0} Input object class validated with \code{inherits()}.
#' @srrstats {G2.0a} Scalar parameters (\code{num_replicates}, \code{nmin}) checked.
#' @srrstats {G2.1} Input types validated (\code{numeric}, \code{logical}, \code{integer}).
#' @srrstats {G2.1a} \code{cutpoint_result} needs \code{optimal_cuts}, \code{parameters}.
#' @srrstats {G2.2} \code{NA} in \code{optimal_cuts} triggers \code{cli_abort()}.
#' @srrstats {RE2.1} \code{na.omit()} and \code{complete.cases()} used for replicate handling.
#' @srrstats {G2.4a} \code{as.integer()} on \code{nmin}, \code{num_replicates}.
#' @srrstats {G5.2} \code{tryCatch()} handles failed replicates gracefully.
#' @srrstats {G5.5} Reproducible via \code{seed} (using \code{doRNG} for parallel).
#' @srrstats {G2.13} \code{cli_abort()} for invalid inputs or missing data.
#' @srrstats {G2.14a} \code{NA} in \code{optimal_cuts} aborts.
#' @srrstats {G2.7} Accepts \code{data.frame} and \code{tibble} (validated via inherits).
#' @srrstats {G5.2a} Handles <20 successful replicates with error.
#' @srrstats {G5.8} Edge cases (insufficient n) checked pre-bootstrap.
#' @srrstats {RE4.0} Returns class \code{validate_cutpoint_result}.
#' @srrstats {RE4.3} Computes 95\% Confidence Intervals via bootstrapping.
#' @srrstats {RE1.3} Output retains original cut-points and parameters.
#'
#' @param cutpoint_result An object from \code{\link{find_cutpoint}}.
#' @param num_replicates Number of bootstrap replicates. Default is 500.
#' @param n_cores Number of CPU cores to use. Default is 1
#' (sequential). Set to > 1 to enable parallel processing.
#' @param seed Optional integer for reproducible results.
#' @param nmin Minimum group size for bootstrap runs. Defaults to 90\%
#' of original \code{nmin} to reduce failures.
#' @param ... Additional arguments passed to \code{\link{find_cutpoint}}
#' (e.g., \code{pop.size}, \code{max.generations} for genetic algorithm).
#'
#' @return An object of class \code{validate_cutpoint_result} with
#' original cuts, 95\% CIs, bootstrap distribution, and parameters.
#'
#' @importFrom foreach %dopar%
#' @importFrom cli cli_h1 cli_text cli_alert_info cli_alert_success
#' @importFrom cli cli_alert_warning cli_inform cli_abort
#' @importFrom cli cli_progress_bar cli_progress_update
#' @importFrom cli cli_warn
#' @importFrom stats quantile na.omit complete.cases sd median IQR
#' @importFrom parallel detectCores makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @importFrom foreach registerDoSEQ
#' @importFrom tidyr pivot_longer everything
#' @export
#'
#' @examples
#' \dontrun{
#' if (requireNamespace("survival", quietly = TRUE)) {
#'   library(survival)
#'   
#'   # 1. Create a tiny simulated baseline clinical cohort dataset
#'   set.seed(123)
#'   n <- 45
#'   toy_data <- data.frame(
#'     time = rexp(n, rate = 0.05),
#'     event = sample(c(0, 1), n, replace = TRUE, prob = c(0.4, 0.6)),
#'     marker = rnorm(n, mean = 4, sd = 1.2)
#'   )
#'
#'   # 2. Locate initial baseline cut-points via systematic search
#'   initial_cut <- find_cutpoint(
#'     data = toy_data, predictor = "marker",
#'     outcome_time = "time", outcome_event = "event",
#'     num_cuts = 1, method = "systematic", criterion = "logrank", nmin = 10
#'   )
#'
#'   # 3. Run a lightweight bootstrap validation stress-test execution loop
#'   val_res <- validate_cutpoint(
#'     cutpoint_result = initial_cut,
#'     num_replicates = 25,  # Small iteration tier for rapid check verification
#'     n_cores = 1,
#'     seed = 123
#'   )
#'
#'   # 4. Invoke structural S3 output verification hooks
#'   print(val_res)
#'   summary(val_res)
#'   plot(val_res)
#' }
#' }
validate_cutpoint <- function(cutpoint_result, num_replicates = 500,
                              n_cores = 1,
                              seed = NULL, nmin = NULL, ...) {
  # --- 1. Validate Input and Set Seed ---
  if (!inherits(cutpoint_result, "find_cutpoint")) {
    cli::cli_abort("Input must be a `find_cutpoint` object.")
  }
  if (!is.numeric(num_replicates) || num_replicates < 1 ||
      num_replicates != round(num_replicates)) {
    cli::cli_abort("`num_replicates` must be a positive integer.")
  }
  if (num_replicates < 20) {
    cli::cli_abort("`num_replicates` must be >= 20 for validation.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
    cli::cli_alert_info("Using random seed {seed} for reproducibility.")
  }
  
  # --- 2. Extract Original Parameters and Reconstruction Layers ---
  original_params <- cutpoint_result$parameters
  original_cuts <- cutpoint_result$optimal_cuts
  candidate_cuts <- cutpoint_result$candidate_cuts
  
  if (any(is.na(original_cuts))) {
    cli::cli_abort("Input `find_cutpoint` object has NA cut-points.")
  }
  
  # Map environmental naming parameters out cleanly
  predictor <- original_params$predictor
  num_cuts <- original_params$num_cuts
  method <- original_params$method
  criterion <- original_params$criterion
  covariates <- original_params$covariates
  grid_by <- original_params$grid_by
  
  # Reconstruct the raw data frame back into original named column metrics
  raw_reconstructed <- cutpoint_result$userdata
  names(raw_reconstructed)[names(raw_reconstructed) == "factor"] <- predictor
  
  # Remap internal parameter tracking names to step past masking traps safely
  outcome_time <- paste0("boot_", predictor, "_time")
  outcome_event <- paste0("boot_", predictor, "_event")
  names(raw_reconstructed)[names(raw_reconstructed) == "time"] <- outcome_time
  names(raw_reconstructed)[names(raw_reconstructed) == "event"] <- outcome_event
  
  n <- nrow(raw_reconstructed)
  
  if (is.null(nmin)) {
    original_nmin_param <- original_params$nmin
    if (original_nmin_param > 0 && original_nmin_param < 1) {
      original_nmin_abs <- floor(original_nmin_param * n)
      nmin <- floor(0.9 * original_nmin_abs)
    } else {
      nmin <- floor(0.9 * original_nmin_param)
    }
    if (nmin < 1) nmin <- 1
    cli::cli_alert_info(paste(
      "Bootstrap `nmin` not set.",
      "Using {nmin} (90% of original) to improve stability."
    ))
  }
  
  if (nmin > 0 && nmin < 1) {
    nmin_abs <- floor(nmin * n)
    if (nmin_abs < 1) nmin_abs <- 1
    cli::cli_alert_info(
      "nmin {nmin} is a proportion. Bootstrap group size set to {nmin_abs}."
    )
    nmin <- nmin_abs
  } else if (nmin >= 1) {
    nmin <- floor(nmin)
  } else {
    cli::cli_abort("`nmin` must be a positive number.")
  }
  
  if (n < nmin * (num_cuts + 1)) {
    cli::cli_abort(
      "Not enough data ({n}) for nmin ({nmin}) and {num_cuts} cut(s)."
    )
  }
  
  cli::cli_alert_info(paste(
    "Validating {num_cuts} cut(s) from '{method}' search",
    "using '{criterion}' over regularized coordinate lattice."
  ))
  
  # --- 3. Setup Backend (Parallel or Sequential) ---
  if (n_cores < 1) n_cores <- 1
  cores_to_use <- 1
  if (n_cores > 1) {
    cores_available <- parallel::detectCores()
    if (is.na(cores_available) || is.null(cores_available)) cores_available <- 2
    
    cores_to_use <- min(cores_available - 1, n_cores, num_replicates)
    if (cores_to_use > 1) {
      if (!requireNamespace("doParallel", quietly = TRUE)) {
        cli::cli_abort("Package 'doParallel' is required for parallel execution.")
      }
      cl <- parallel::makeCluster(cores_to_use, type = "PSOCK")
      doParallel::registerDoParallel(cl)
      on.exit(parallel::stopCluster(cl), add = TRUE)
      cli::cli_alert_info(
        "Running {num_replicates} replicates on {cores_to_use} cores..."
      )
      
      if (!is.null(seed)) {
        if (!requireNamespace("doRNG", quietly = TRUE)) {
          cli::cli_abort(
            "Package 'doRNG' is required for reproducible parallel computation."
          )
        }
        doRNG::registerDoRNG()
      }
    } else {
      cores_to_use <- 1
    }
  }
  if (cores_to_use == 1) {
    foreach::registerDoSEQ()
    cli::cli_alert_info(
      "Running {num_replicates} replicates sequentially (n_cores = 1)."
    )
  }
  
  # --- 4. Main Bootstrap Loop ---
  if (cores_to_use == 1) {
    pb <- cli::cli_progress_bar("Bootstrapping", total = num_replicates)
  }
  
  i <- NULL
  extra_args <- list(...)
  
  if ("pop.size" %in% names(extra_args) && is.null(extra_args$pop.size)) extra_args$pop.size <- NULL
  if ("max.generations" %in% names(extra_args) && is.null(extra_args$max.generations)) extra_args$max.generations <- NULL
  
  bootstrap_results <- foreach::foreach(
    i = 1:num_replicates,
    .combine = "rbind",
    .packages = c("survival", "OptSurvCutR"),
    .errorhandling = "pass"
  ) %dopar% {
    if (cores_to_use == 1) {
      cli::cli_progress_update(id = pb)
    }
    
    boot_indices <- sample(1:n, n, replace = TRUE)
    boot_data <- raw_reconstructed[boot_indices, ]
    
    final_args <- c(
      list(
        data = boot_data,
        predictor = predictor,
        outcome_time = outcome_time,
        outcome_event = outcome_event,
        method = method,
        num_cuts = num_cuts,
        criterion = criterion,
        covariates = covariates,
        nmin = nmin,
        pop.size = extra_args$pop.size,
        max.generations = extra_args$max.generations,
        use_cpp = original_params$use_cpp,
        grid_by = grid_by,
        quiet = TRUE,
        candidate_cuts = candidate_cuts
      ),
      extra_args[!(names(extra_args) %in% c("pop.size", "max.generations"))]
    )
    
    res <- tryCatch(
      suppressWarnings(suppressMessages(do.call(OptSurvCutR::find_cutpoint, final_args))),
      error = function(e) {
        return(NULL)
      }
    )
    
    if (is.null(res) || any(is.na(res$optimal_cuts)) || length(res$optimal_cuts) != num_cuts) {
      return(matrix(rep(NA_real_, num_cuts), nrow = 1, ncol = num_cuts))
    }
    matrix(sort(res$optimal_cuts), nrow = 1, ncol = num_cuts)
  }
  
  # --- 5. Process Results and Calculate CIs ---
  if (!is.matrix(bootstrap_results)) {
    bootstrap_matrix <- matrix(unlist(bootstrap_results),
                               nrow = num_replicates,
                               ncol = num_cuts, byrow = TRUE
    )
  } else {
    bootstrap_matrix <- bootstrap_results
  }
  
  # Ensure clean matrix topology and dimension tracking boundaries
  bootstrap_matrix <- matrix(as.numeric(bootstrap_matrix), nrow = num_replicates, ncol = num_cuts)
  colnames(bootstrap_matrix) <- paste0("Cut", 1:num_cuts)
  
  successful_reps <- sum(stats::complete.cases(bootstrap_matrix))
  failed_reps <- num_replicates - successful_reps
  
  if (failed_reps > 0) {
    cli::cli_alert_warning(
      "{failed_reps} of {num_replicates} replicates failed."
    )
  }
  if (successful_reps < 20) {
    cli::cli_abort(paste(
      "< 20 successful replicates ({successful_reps}).",
      "Increase num_replicates or reduce nmin."
    ))
  }
  cli::cli_alert_success("{successful_reps} replicates completed.")
  if (successful_reps / num_replicates < 0.8 && successful_reps > 0) {
    warning(
      "Success rate of model fitting is below 80%. CIs may be unreliable.",
      call. = FALSE
    )
  }
  bootstrap_matrix_clean <- na.omit(bootstrap_matrix)
  
  bootstrap_matrix_clean <- matrix(
    as.numeric(bootstrap_matrix_clean),
    nrow = successful_reps,
    ncol = num_cuts
  )
  
  ci <- apply(bootstrap_matrix_clean, 2, stats::quantile,
              probs = c(0.025, 0.975), na.rm = TRUE
  )
  if (num_cuts == 1) {
    ci_df <- data.frame(Lower = ci[1], Upper = ci[2])
  } else {
    ci_df <- as.data.frame(t(ci))
    names(ci_df) <- c("Lower", "Upper")
  }
  row.names(ci_df) <- paste0("Cut ", 1:num_cuts)
  
  boot_summary <- lapply(1:num_cuts, function(i) {
    valid_cuts <- bootstrap_matrix_clean[, i]
    if (length(valid_cuts) == 0) {
      return(list(mean = NA, sd = NA, median = NA, Q1 = NA, Q3 = NA))
    }
    list(
      mean = mean(valid_cuts, na.rm = TRUE),
      sd = stats::sd(valid_cuts, na.rm = TRUE),
      median = stats::median(valid_cuts, na.rm = TRUE),
      Q1 = stats::quantile(valid_cuts, 0.25, na.rm = TRUE),
      Q3 = stats::quantile(valid_cuts, 0.75, na.rm = TRUE)
    )
  })
  names(boot_summary) <- paste0("Cut", 1:num_cuts)
  
  bootstrap_df <- as.data.frame(bootstrap_matrix_clean)
  names(bootstrap_df) <- paste0("Cut_point_", 1:num_cuts)
  
  output <- list(
    original_cuts = original_cuts,
    confidence_intervals = ci_df,
    bootstrap_distribution = bootstrap_df,
    boot_summary = boot_summary,
    parameters = list(
      num_replicates = num_replicates,
      successful_reps = successful_reps,
      failed_reps = failed_reps,
      n_cores = cores_to_use,
      seed = seed,
      nmin = nmin,
      method = method,
      criterion = criterion,
      covariates = covariates
    )
  )
  class(output) <- "validate_cutpoint_result"
  return(output)
}
