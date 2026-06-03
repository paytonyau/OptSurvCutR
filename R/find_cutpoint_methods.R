# ===================================================================
# S3 METHODS
# Print, Summary, and Plot methods for all 'find_cutpoint' objects
# ===================================================================

#' @rdname find_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.17} `print()` method provided.
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
                       "Optimal Statistic") # Fallback

  if (x$parameters$criterion == "p_value" && x$parameters$method == "genetic") {
    stat_label <- "Optimal LRT Statistic"
  }

  # Safely handle missing optimal stats
  stat_val_fmt <- if (!is.null(x$optimal_stat) && is.numeric(x$optimal_stat)) {
    round(x$optimal_stat, 4)
  } else {
    "N/A"
  }

  rounded_cuts <- round(x$optimal_cuts, 3)

  bullets <- c(
    "*" = "Predictor: {.strong {x$parameters$predictor}}",
    "*" = "Criterion: {.strong {x$parameters$criterion}}",
    "*" = "{stat_label}: {.strong {stat_val_fmt}}",
    "v" = "Recommended Cut-point(s): {.strong {paste(rounded_cuts, collapse = ', ')}}"
  )

  if (!is.null(x$permuted_p_value) && !is.na(x$permuted_p_value)) {
    bullets <- c(bullets, "*" = "Permuted P-value ({x$n_perm} runs): {.strong {round(x$permuted_p_value, 4)}}")
  }

  cli::cli_bullets(bullets)
  cli::cli_text("\nHint: Use `summary()` for clinical details and Cox regression.")
  invisible(x)
}

#' @rdname find_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.18} `summary()` method provided.
#' @srrstats {RE1.1} Assumes PH; check `summary()` for `cox.zph`.
#' @srrstats {RE1.4} Cox PH assumption test via `summary(fit)$cox_zph`.
#' @srrstats {RE2.1} Estimates and standard errors available from `coxph` in `summary()`.
#' @export
summary.find_cutpoint <- function(object, show_model = TRUE, show_group_counts = TRUE,
                                  show_medians = TRUE, show_ph_test = TRUE, show_params = TRUE, ...) {

  # --- HEADER & CONCLUSION ---
  cli::cli_h1(paste("Optimal Cut-point Analysis for Survival Data ({tools::toTitleCase(object$parameters$method)})"))

  # --- NEW: THE TROUBLESHOOTING FAILURE CATCHER ---
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
      cli::cli_bullets(c(
        ">" = "Please double-check your settings.",
        ">" = "{.strong Relax constraints:} Lower {.code nmin} (e.g., 0.05) or set {.code boundary.enforcement = 1}.",
        ">" = "{.strong Check feasibility:} Run with {.code method = 'systematic'} to verify a solution physically exists.",
        ">" = "{.strong Reduce complexity:} Try searching for {.val {max(1, object$parameters$num_cuts - 1)}} cut-point(s) instead."
      ))
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

  # CRAN SAFEGUARD: Ensure cuts are unique to prevent `cut()` breaking
  safe_cuts <- unique(object$optimal_cuts)
  num_safe_groups <- length(safe_cuts) + 1

  data$group <- factor(cut(data$factor, breaks = c(-Inf, safe_cuts, Inf), labels = paste0("G", 1:num_safe_groups)))

  # --- 1. CLINICAL RISK COHORTS ---
  if (show_group_counts || show_medians) {
    cli::cli_h2("1. Stratified Risk Cohorts")

    counts_table <- as.data.frame(table(data$group))
    names(counts_table) <- c("Group", "N")

    # Safely aggregate events (prevents groups with 0 events from vanishing)
    event_counts <- stats::aggregate(event ~ group, data = data, FUN = sum, drop = FALSE, na.action = stats::na.pass)
    names(event_counts) <- c("Group", "Events")

    cohort_df <- merge(counts_table, event_counts, by = "Group", all.x = TRUE)
    cohort_df$Events[is.na(cohort_df$Events)] <- 0

    if (show_medians) {
      fit_km <- survival::survfit(survival::Surv(time, event) ~ group, data = data)
      surv_summary <- summary(fit_km)$table

      # CRAN SAFEGUARD: Robust Merge mapping (prevents matrix mismatch if survfit drops empty cohorts)
      if (is.null(dim(surv_summary))) {
        surv_df <- as.data.frame(t(surv_summary))
        # Find the first group that actually has data
        surv_df$Group <- as.character(cohort_df$Group[cohort_df$N > 0][1])
      } else {
        surv_df <- as.data.frame(surv_summary)
        surv_df$Group <- gsub("group=", "", rownames(surv_df))
      }

      # Merge safely by Group Name
      cohort_df <- merge(cohort_df, surv_df[, c("Group", "median", "0.95LCL", "0.95UCL"), drop = FALSE], by = "Group", all.x = TRUE)

      cohort_df$Median_CI <- sprintf(
        "%s (%s - %s)",
        ifelse(is.na(cohort_df$median), "NA", round(cohort_df$median, 1)),
        ifelse(is.na(cohort_df$`0.95LCL`), "NA", round(cohort_df$`0.95LCL`, 1)),
        ifelse(is.na(cohort_df$`0.95UCL`), "NA", round(cohort_df$`0.95UCL`, 1))
      )

      print(cohort_df[, c("Group", "N", "Events", "Median_CI")], row.names = FALSE, right = TRUE)
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

  # --- 2. COX PROPORTIONAL-HAZARDS ---
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
          Group = gsub("group", "", rownames(coefs)),
          HR = round(conf[, "exp(coef)"], 3),
          Lower = round(conf[, "lower .95"], 3),
          Upper = round(conf[, "upper .95"], 3),
          P_Value = round(coefs[, "Pr(>|z|)"], 3)
        )

        cox_df$Signif <- as.character(symnum(cox_df$P_Value, corr = FALSE, na = FALSE,
                                             cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                                             symbols = c("***", "**", "*", ".", " ")))

        print(cox_df, row.names = FALSE, right = TRUE)
        cat("\n")
      }

      conc <- round(cox_sum$concordance["C"], 3)
      pval <- round(cox_sum$sctest["pvalue"], 3)
      cli::cli_alert_info("Overall Model: Concordance = {.val {conc}} | Log-rank p = {.val {pval}}")
      cat("\n")
    }
  }

  # --- 3. TIME-DEPENDENT DIAGNOSTICS (Schoenfeld) ---
  if (show_ph_test && !is.null(fit_cox)) {
    cli::cli_h2("3. Time-Dependent Diagnostics (Schoenfeld)")

    zph_test <- tryCatch(survival::cox.zph(fit_cox), error = function(e) NULL)

    if (!is.null(zph_test)) {
      zph_table <- zph_test$table

      # SAFE EXTRACTION: Handle 1-cut vs multi-cut GLOBAL row differences
      global_p <- if ("GLOBAL" %in% rownames(zph_table)) {
        zph_table["GLOBAL", "p"]
      } else {
        zph_table[nrow(zph_table), "p"] # Fallback for 1-cut models
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

  # --- 4. PARAMETERS ---
  if (show_params) {
    cli::cli_h2("4. Analysis Parameters")
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
    if (!is.null(object$n_perm) && object$n_perm > 0) {
      param_bullets <- c(param_bullets, "*" = "Permutations: {object$n_perm}")
    }
    cli::cli_bullets(param_bullets)
    cat("\n")
  }

  invisible(object)
}

#' @rdname find_cutpoint
#' @section srrstats compliance:
#' .
#' @srrstats {RE6.0} Primary plot method for the package class.
#' @srrstats {RE6.1} Generates standard survival curves (outcome).
#' @srrstats {G2.4b} Strict `match.arg` used for plot types.
#' @srrstats {G5.3} Provides an escape hatch via `return_data = TRUE`.
#' @export
plot.find_cutpoint <- function(x, type = c("outcome", "distribution", "forest", "diagnostic", "auc", "all"),
                               return_data = FALSE, ...) {
  type <- match.arg(type)

  if (any(is.na(x$optimal_cuts))) {
    message("Cannot generate plot: no valid cut-point")
    return(invisible(NULL))
  }

  df <- x$userdata
  cuts <- sort(x$optimal_cuts)
  labels <- paste0("G", 1:(length(cuts) + 1))
  df$group <- factor(findInterval(df$factor, cuts, left.open = TRUE) + 1L, labels = labels)

  args <- list(...)
  valid_ref <- FALSE
  ref_name <- ""

  if ("reference_group" %in% names(args)) {
    ref <- args$reference_group
    if (!is.null(ref) && ref %in% levels(df$group)) {
      df$group <- stats::relevel(df$group, ref = ref)
      valid_ref <- TRUE
      ref_name <- ref
    } else {
      message("Invalid reference group. Defaulting to baseline group.")
    }
    args$reference_group <- NULL
  }

  if (return_data) return(df)

  if (type == "diagnostic") return(plot_cutpoint_residuals(x, ...))
  if (type == "auc") return(plot_time_dependent_auc(x))

  km_form_str <- "survival::Surv(time, event) ~ group"
  km_form <- stats::as.formula(km_form_str)

  cox_form_str <- km_form_str
  if (!is.null(x$parameters$covariates)) {
    cox_form_str <- paste(cox_form_str, "+", paste(x$parameters$covariates, collapse = " + "))
  }
  cox_form <- stats::as.formula(cox_form_str)

  if (type == "distribution") {
    p_dist <- ggplot2::ggplot(df, ggplot2::aes(x = .data$factor)) +
      ggplot2::geom_density(fill = "#0072B2", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = cuts, color = "#D55E00", linetype = "dashed", linewidth = 1) +
      ggplot2::labs(title = "Predictor Distribution", subtitle = paste("Cut-points:", paste(round(cuts, 3), collapse = ", ")),
                    x = "Predictor Value", y = "Density") +
      theme_optsurv()
    return(p_dist)
  }

  if (type == "outcome") {
    fit <- survival::survfit(km_form, data = df)
    fit$call$formula <- parse(text = km_form_str)[[1]]

    if (is.null(args$palette)) args$palette <- "nejm"
    if (is.null(args$pval)) args$pval <- TRUE
    if (is.null(args$risk.table)) args$risk.table <- TRUE

    surv_args <- c(list(fit = fit, data = df, ggtheme = theme_optsurv()), args)
    return(do.call(survminer::ggsurvplot, surv_args))
  }

  if (type == "forest") {
    fit_cox <- tryCatch(survival::coxph(cox_form, data = df), error = function(e) NULL)
    if (is.null(fit_cox)) {
      message("Could not fit Cox model for forest plot")
      return(invisible(NULL))
    }
    fit_cox$call$formula <- parse(text = cox_form_str)[[1]]

    forest_args <- c(list(model = fit_cox, data = df), args)
    p_forest <- do.call(survminer::ggforest, forest_args)

    if (valid_ref) {
      p_forest <- p_forest + ggplot2::labs(x = "Hazard Ratio (95% CI)", subtitle = paste("Reference group:", ref_name))
    } else {
      p_forest <- p_forest + ggplot2::labs(x = "Hazard Ratio (95% CI)")
    }
    return(p_forest)
  }

  if (type == "all") {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      cli::cli_abort(c("Package {.pkg patchwork} is required for the dashboard view.", "i" = "Install it using: install.packages('patchwork')"))
    }
    p_dist <- ggplot2::ggplot(df, ggplot2::aes(x = .data$factor)) +
      ggplot2::geom_density(fill = "#0072B2", alpha = 0.5) +
      ggplot2::geom_vline(xintercept = cuts, color = "#D55E00", linetype = "dashed", linewidth = 1) +
      ggplot2::labs(title = "Predictor Distribution", subtitle = paste("Cut-points:", paste(round(cuts, 3), collapse = ", ")),
                    x = "Predictor Value", y = "Density") +
      theme_optsurv()

    fit <- survival::survfit(km_form, data = df)
    fit$call$formula <- parse(text = km_form_str)[[1]]
    if (is.null(args$palette)) args$palette <- "nejm"
    if (is.null(args$pval)) args$pval <- TRUE
    args$risk.table <- FALSE

    surv_args <- c(list(fit = fit, data = df, ggtheme = theme_optsurv()), args)
    p_surv <- do.call(survminer::ggsurvplot, surv_args)

    dashboard <- p_dist / p_surv$plot + patchwork::plot_layout(heights = c(1, 2))
    return(dashboard)
  }
}
