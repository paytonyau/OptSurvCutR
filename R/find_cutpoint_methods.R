# ===================================================================
# S3 METHODS
# Print and Summary methods for all 'find_cutpoint' objects
# ===================================================================

#' @rdname find_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.17} \code{print()} method provided.
#' @param x A find_cutpoint result object.
#' @param ... Additional arguments passed down to downstream rendering pipelines.
#' @export
print.find_cutpoint <- function(x, ...) {
  if (is.null(x)) {
    cli::cli_inform("No cut-point object.")
    return(invisible(x))
  }
  
  if (all(is.na(x$optimal_cuts))) {
    cli::cli_inform("No optimal cut-point determined.")
    return(invisible(x))
  }
  
  method_name <- tools::toTitleCase(x$parameters$method)
  cli::cli_h1(paste("Optimal Cut-point Analysis for Survival Data ({method_name})"))
  
  stat_label <- switch(x$parameters$criterion,
                       "logrank" = "Optimal Log-Rank Statistic",
                       "hazard_ratio" = "Optimal Hazard Ratio",
                       "p_value" = "Optimal P-value",
                       "Optimal Statistic"
  )
  
  if (x$parameters$criterion == "p_value" && x$parameters$method == "genetic") {
    stat_label <- "Optimal LRT Statistic"
  }
  
  stat_val_fmt <- if (!is.null(x$optimal_stat) && is.numeric(x$optimal_stat)) {
    round(x$optimal_stat, 4)
  } else {
    "N/A"
  }
  
  rounded_cuts <- round(x$optimal_cuts, 3)
  
  # CRAN/goodpractice FIX: Create an unnamed vector first, then set names explicitly 
  # to prevent duplicated argument errors during the c() collection parse step.
  bullets <- c(
    paste0("Predictor: {.strong ", x$parameters$predictor, "}"),
    paste0("Criterion: {.strong ", x$parameters$criterion, "}"),
    paste0(stat_label, ": {.strong ", stat_val_fmt, "}"),
    paste0("Recommended Cut-point(s): {.strong ", paste(rounded_cuts, collapse = ', '), "}")
  )
  names(bullets) <- c("*", "*", "*", "v")
  
  if (!is.null(x$permuted_p_value) && !is.na(x$permuted_p_value)) {
    new_bullet <- paste0("Permuted P-value (", x$n_perm, " runs): {.strong ", round(x$permuted_p_value, 4), "}")
    names(new_bullet) <- "*"
    bullets <- c(bullets, new_bullet)
  }
  
  cli::cli_bullets(bullets)
  cli::cli_text("\nHint: Use `summary()` for clinical details and Cox regression.")
  invisible(x)
}

#' @rdname find_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.18} \code{summary()} method provided.
#' @srrstats {RE1.1} Assumes PH; check \code{summary()} for \code{cox.zph}.
#' @srrstats {RE1.4} Cox PH assumption test via \code{summary(fit)$cox_zph}.
#' @srrstats {RE2.1} Estimates and standard errors available from \code{coxph} in \code{summary()}.
#' @param object A find_cutpoint result object for summary evaluation.
#' @param show_model Logical. Whether to print the full Cox model summary frame.
#' @param show_group_counts Logical. Whether to show stratified sample split counts.
#' @param show_medians Logical. Whether to display Kaplan-Meier median tracking times.
#' @param show_ph_test Logical. Display the proportional hazards validation check.
#' @param show_params Logical. Print original baseline parameters.
#' @param ... Additional arguments passed down to downstream rendering pipelines.
#' @importFrom stats as.formula aggregate symnum time na.pass
#' @importFrom survival Surv coxph cox.zph survfit
#' @export
summary.find_cutpoint <- function(object, show_model = TRUE, show_group_counts = TRUE,
                                  show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE, ...) {
  cli::cli_h1(paste("Optimal Cut-point Analysis for Survival Data ({tools::toTitleCase(object$parameters$method)})"))
  
  if (is.null(object) || is.null(object$optimal_cuts) || all(is.na(object$optimal_cuts))) {
    cli::cli_alert_danger("No valid optimal cut-point found.")
    cat("\n")
    
    if (!is.null(object$parameters)) {
      cli::cli_alert_info(
        "The algorithm failed to find a solution for {.val {object$parameters$num_cuts}} cut-points. ",
        "This occurs when your search constraints conflict with the data distribution."
      )
      cat("\n")
      cli::cli_alert_success("{.strong Recommended fixes:}")
      
      fix_bullets <- c(
        "Please double-check your settings.",
        "{.strong Relax constraints:} Lower {.code nmin} (e.g., 0.05) or set {.code boundary.enforcement = 1}.",
        "{.strong Check feasibility:} Run with {.code method = 'systematic'} to verify a solution physically exists.",
        paste0("{.strong Reduce complexity:} Try searching for {.val ", max(1, object$parameters$num_cuts - 1), "} cut-point(s) instead.")
      )
      names(fix_bullets) <- c(">", ">", ">", ">")
      cli::cli_bullets(fix_bullets)
    }
    return(invisible(object))
  }
  
  best_cuts <- round(object$optimal_cuts, 3)
  cli::cli_alert_success("Optimal Threshold(s): {.val {paste(best_cuts, collapse = ', ')}}")
  
  if (!is.null(object$permuted_p_value) && !is.na(object$permuted_p_value)) {
    cli::cli_alert_info("Permutation-Adjusted P-value ({object$n_perm} runs): {.val {round(object$permuted_p_value, 4)}}")
  }
  cat("\n")
  
  data <- object$userdata
  
  safe_cuts <- unique(object$optimal_cuts)
  num_safe_groups <- length(safe_cuts) + 1
  
  data$group <- factor(cut(data$factor, breaks = c(-Inf, safe_cuts, Inf), labels = paste0("G", 1:num_safe_groups)))
  
  if (show_group_counts || show_medians) {
    cli::cli_h2("1. Stratified Risk Cohorts")
    
    counts_table <- as.data.frame(table(data$group))
    names(counts_table) <- c("Group", "N")
    
    event_counts <- stats::aggregate(event ~ group, data = data, FUN = sum, drop = FALSE, na.action = stats::na.pass)
    names(event_counts) <- c("Group", "Events")
    
    cohort_df <- merge(counts_table, event_counts, by = "Group", all.x = TRUE)
    cohort_df$Events[is.na(cohort_df$Events)] <- 0
    
    if (show_medians) {
      fit_km <- survival::survfit(survival::Surv(time, event) ~ group, data = data)
      surv_summary <- summary(fit_km)$table
      
      if (is.null(dim(surv_summary))) {
        surv_df <- as.data.frame(t(surv_summary))
        surv_df$Group <- as.character(cohort_df$Group[cohort_df$N > 0][1])
      } else {
        surv_df <- as.data.frame(surv_summary)
        surv_df$Group <- gsub("group=", "", rownames(surv_df), fixed = TRUE)
      }
      
      horizon_time <- max(data$time, na.rm = TRUE)
      summary_horizon <- summary(fit_km, times = horizon_time, extend = TRUE)
      
      landmark_df <- data.frame(
        Group = gsub("group=", "", summary_horizon$strata, fixed = TRUE),
        Surv_Rate = paste0(
          round(summary_horizon$surv * 100, 1), "% (",
          round(summary_horizon$lower * 100, 1), "-",
          round(summary_horizon$upper * 100, 1), ")"
        )
      )
      
      cohort_df <- merge(cohort_df, surv_df[, c("Group", "median"), drop = FALSE], by = "Group", all.x = TRUE)
      cohort_df <- merge(cohort_df, landmark_df, by = "Group", all.x = TRUE)
      
      cohort_df$Median_Time <- ifelse(is.na(cohort_df$median), "NR (Not Reached)", as.character(round(cohort_df$median, 1)))
      names(cohort_df)[names(cohort_df) == "Surv_Rate"] <- paste0("OS Rate at T=", round(horizon_time, 1))
      
      print(cohort_df[, c("Group", "N", "Events", "Median_Time", paste0("OS Rate at T=", round(horizon_time, 1)))], row.names = FALSE, right = TRUE)
    } else {
      print(cohort_df, row.names = FALSE, right = TRUE)
    }
    cat("\n")
  }
  
  formula_str <- "survival::Surv(time, event) ~ group"
  if (!is.null(object$parameters$covariates)) {
    formula_str <- paste(formula_str, "+", paste(object$parameters$covariates, collapse = " + "))
  }
  fit_cox <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = data), error = function(e) NULL)
  
  if (show_model) {
    cli::cli_h2("2. Cox Proportional-Hazards")
    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model: convergence failed.")
    } else {
      cox_sum <- summary(fit_cox)
      coefs <- cox_sum$coefficients
      conf <- cox_sum$conf.int
      
      if (nrow(coefs) > 0) {
        cox_df <- data.frame(
          Group = gsub("group", "", rownames(coefs), fixed = TRUE),
          HR = round(conf[, "exp(coef)"], 3),
          Lower = round(conf[, "lower .95"], 3),
          Upper = round(conf[, "upper .95"], 3),
          P_Value = round(coefs[, "Pr(>|z|)"], 3)
        )
        
        cox_df$Signif <- as.character(stats::symnum(cox_df$P_Value,
                                                    corr = FALSE, na = FALSE,
                                                    cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                                    symbols = c("***", "**", "*", ".", " ")
        ))
        
        print(cox_df, row.names = FALSE, right = TRUE)
        cat("\n")
      }
      
      conc <- round(cox_sum$concordance["C"], 3)
      pval <- round(cox_sum$sctest["pvalue"], 3)
      cli::cli_alert_info("Overall Model: Concordance = {.val {conc}} | Log-rank p = {.val {pval}}")
      cat("\n")
    }
  }
  
  if (show_ph_test && !is.null(fit_cox)) {
    cli::cli_h2("3. Time-Dependent Diagnostics (Schoenfeld)")
    
    zph_test <- tryCatch(survival::cox.zph(fit_cox), error = function(e) NULL)
    
    if (!is.null(zph_test)) {
      zph_table <- zph_test$table
      
      global_p <- if ("GLOBAL" %in% rownames(zph_table)) {
        zph_table["GLOBAL", "p"]
      } else {
        zph_table[nrow(zph_table), "p"]
      }
      
      if (!is.na(global_p) && global_p < 0.05) {
        cli::cli_alert_info("Insight: The hazard ratios appear to fluctuate over time (Global p = {.val {signif(global_p, 3)}}).")
        cli::cli_text("The cut-points successfully separate the subjects, but the relative event risk between these cohorts likely evolves as follow-up time increases.")
      } else if (!is.na(global_p)) {
        cli::cli_alert_success("Passed: The proportional hazards assumption holds across the follow-up period (Global p = {.val {signif(global_p, 3)}}).")
      } else {
        cli::cli_text("Could not compute Global PH assumption test.")
      }
    } else {
      cli::cli_text("Could not compute PH assumption test.")
    }
    cat("\n")
  }
  
  if (show_params) {
    cli::cli_h2("4. Analysis Parameters")
    params <- object$parameters
    
    param_bullets <- c(
      paste0("Search Method: ", tools::toTitleCase(params$method)),
      paste0("Predictor: ", params$predictor),
      paste0("Number of cuts: ", params$num_cuts),
      paste0("Minimum group size (nmin): ", params$nmin)
    )
    names(param_bullets) <- c("*", "*", "*", "*")
    
    if (!is.null(params$covariates)) {
      cov_bullet <- paste0("Covariates: ", paste(params$covariates, collapse = ', '))
      names(cov_bullet) <- "*"
      param_bullets <- c(param_bullets, cov_bullet)
    }
    if (!is.null(object$n_perm) && object$n_perm > 0) {
      perm_bullet <- paste0("Permutations: ", object$n_perm)
      names(perm_bullet) <- "*"
      param_bullets <- c(param_bullets, perm_bullet)
    }
    cli::cli_bullets(param_bullets)
    cat("\n")
  }
  
  invisible(object)
}
