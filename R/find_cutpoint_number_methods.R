# ===================================================================
# S3 METHODS
# Print, Summary, and Plot methods for 'find_cutpoint_number_result'
# ===================================================================

#' @rdname find_cutpoint_number
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.17} `print()` method provided.
#' @export
print.find_cutpoint_number_result <- function(x, ...) {
  cli::cli_h1("Optimal Cut-point Number Analysis")

  method_text <- x$parameters$method %||% "Unknown"
  criterion_text <- x$parameters$criterion %||% "IC"

  cli::cli_text("Method: {.strong {tools::toTitleCase(method_text)}}")
  cli::cli_text("Criterion: {.strong {toupper(criterion_text)}}")

  if (!is.null(x$parameters$covariates)) {
    cli::cli_text("Covariates: {.strong {paste(x$parameters$covariates, collapse = ', ')}}")
  }

  if (is.null(x$results) || nrow(x$results) == 0 || all(is.na(x$results[[criterion_text]]))) {
    cli::cli_inform("No optimal model could be determined.")
    return(invisible(x))
  }

  print_df <- x$results

  # Clean the cuts column list into readable strings
  if ("cuts" %in% names(print_df) && is.list(print_df$cuts)) {
    print_df$cuts <- vapply(print_df$cuts, function(c) {
      if (is.null(c)) "NA" else paste(round(c, 3), collapse = ", ")
    }, FUN.VALUE = character(1))
  }

  # Round numeric columns safely
  is_num <- vapply(print_df, is.numeric, FUN.VALUE = logical(1))
  print_df[is_num] <- lapply(print_df[is_num], round, 2)

  # Format Weights to percentages safely (bypassing the rounding above)
  weight_col <- names(print_df)[grepl("_Weight$", names(print_df))]
  if (length(weight_col) > 0 && weight_col[1] %in% names(x$results)) {
    target_col <- weight_col[1]
    print_df[[target_col]] <- paste0(round(x$results[[target_col]] * 100, 1), "%")
  }

  # Extract and print table
  final_cols <- c("num_cuts", criterion_text, paste0("Delta_", criterion_text), weight_col[1], "Evidence", "cuts")
  final_cols_exist <- final_cols[final_cols %in% names(print_df)]
  print(print_df[, final_cols_exist, drop = FALSE], row.names = FALSE, right = TRUE)

  # Find and announce the winner
  best_result <- x$results[which.min(x$results[[criterion_text]]), ]

  if (nrow(best_result) > 0 && is.finite(best_result[[criterion_text]])) {
    best_cuts_vals <- best_result$cuts[[1]]
    cli::cli_alert_success(paste(
      "\nConclusion: {best_result$num_cuts} cut-point(s) is",
      "best based on {toupper(criterion_text)}."
    ))
    if (!is.null(best_cuts_vals)) {
      rounded_cuts <- round(best_cuts_vals, 3)
      cli::cli_text("  Optimal cuts at: {.strong {paste(rounded_cuts, collapse = ', ')}}")
    }
  } else {
    cli::cli_inform("\nConclusion: No optimal model could be determined.")
  }

  cli::cli_text("\nHint: Use `summary()` for clinical details, PH diagnostics, and Cox regression.")
  invisible(x)
}


#' @param show_comparison_table Logical. Show model comparison table?
#' @param show_best_model_details Logical. Show details for best model?
#' @param show_group_counts Logical. Show group counts for best model?
#' @param show_medians Logical. Show median survival for best model?
#' @param show_ph_test Logical. Show Proportional Hazards assumption test?
#' @param plot.it Logical. Display model selection plot?
#' @rdname find_cutpoint_number
#' @section srrstats compliance:
#' .
#' @srrstats {RE4.18} `summary()` method provided.
#' @srrstats {RE1.3} PH diagnostics via `cox.zph` in `summary()`.
#' @srrstats {RE1.3a} `summary()` includes PH test results.
#' @srrstats {RE2.0} Estimates/SEs from `coxph` in `summary()`.
#' @export
summary.find_cutpoint_number_result <- function(
    object, show_comparison_table = TRUE,
    show_best_model_details = TRUE,
    show_group_counts = TRUE, show_medians = TRUE,
    show_ph_test = TRUE, plot.it = FALSE, ...
) {
  criterion_text <- object$parameters$criterion %||% "IC"

  cli::cli_h1(paste("Optimal Cut-point Number Analysis ({tools::toTitleCase(object$parameters$method %||% \"Unknown\")})"))

  if (is.null(object) || is.null(object$results) ||
      nrow(object$results) == 0 ||
      all(is.na(object$results[[criterion_text]]))) {
    cli::cli_inform("Cannot summarise: no valid model was found.")
    return(invisible(object))
  }

  best_result <- object$results[which.min(object$results[[criterion_text]]), ]
  num_cuts <- best_result$num_cuts
  best_cuts_vals <- best_result$cuts[[1]]

  # --- HEADER & CONCLUSION ---
  cli::cli_alert_success("Best Model: {.strong {num_cuts} Cut-points} (Criterion: {toupper(criterion_text)})")
  if (!is.null(best_cuts_vals)) {
    cli::cli_alert_info("Optimal Thresholds: {.val {paste(round(best_cuts_vals, 3), collapse = ', ')}}")
  }
  cat("\n")

  # --- 1. MODEL COMPARISON TABLE (Deduplicated) ---
  if (show_comparison_table) {
    cli::cli_h2("1. Model Comparison")

    comp_df <- object$results
    comp_df$cuts <- NULL # Remove list column for printing

    # Add winner arrow
    comp_df$Marker <- ifelse(comp_df$num_cuts == num_cuts, ">", " ")

    weight_col <- names(comp_df)[grepl("_Weight$", names(comp_df))]
    if (length(weight_col) > 0) {
      target_col <- weight_col[1]
      comp_df[[target_col]] <- paste0(round(comp_df[[target_col]] * 100, 1), "%")
    }

    cols_to_print <- c("Marker", "num_cuts", criterion_text, paste0("Delta_", criterion_text), weight_col[1], "Evidence")
    cols_exist <- cols_to_print[cols_to_print %in% names(comp_df)]

    is_num <- vapply(comp_df[, cols_exist, drop=FALSE], is.numeric, FUN.VALUE = logical(1))
    comp_df[, cols_exist][is_num] <- lapply(comp_df[, cols_exist][is_num], round, 2)

    print(comp_df[, cols_exist, drop = FALSE], row.names = FALSE, right = TRUE)
    cat("\n")
  }

  # --- 2. CLINICAL RISK COHORTS ---
  if (show_best_model_details) {
    data <- object$userdata
    cov_part <- if (!is.null(object$parameters$covariates)) paste(" +", paste(object$parameters$covariates, collapse = " + ")) else ""

    if (num_cuts > 0) {
      # CRAN SAFEGUARD: Unique cuts to prevent breaks error
      safe_cuts <- unique(best_cuts_vals)
      num_safe_groups <- length(safe_cuts) + 1

      data$group <- factor(cut(data$factor, breaks = c(-Inf, safe_cuts, Inf), labels = paste0("G", 1:num_safe_groups)))

      if (show_group_counts || show_medians) {
        cli::cli_h2("2. Clinical Risk Cohorts")

        counts_table <- as.data.frame(table(data$group))
        names(counts_table) <- c("Group", "N")

        # Zero-event safeguard
        event_counts <- stats::aggregate(event ~ group, data = data, FUN = sum, drop = FALSE, na.action = stats::na.pass)
        names(event_counts) <- c("Group", "Events")

        cohort_df <- merge(counts_table, event_counts, by = "Group", all.x = TRUE)
        cohort_df$Events[is.na(cohort_df$Events)] <- 0

        if (show_medians) {
          fit_km <- survival::survfit(survival::Surv(time, event) ~ group, data = data)
          surv_summary <- summary(fit_km)$table

          # Matrix drop safeguard
          if (is.null(dim(surv_summary))) {
            surv_df <- as.data.frame(t(surv_summary))
            surv_df$Group <- as.character(cohort_df$Group[cohort_df$N > 0][1])
          } else {
            surv_df <- as.data.frame(surv_summary)
            surv_df$Group <- gsub("group=", "", rownames(surv_df))
          }

          # CRAN SAFEGUARD: Only extract columns that successfully generated
          target_cols <- intersect(c("Group", "median", "0.95LCL", "0.95UCL"), names(surv_df))
          cohort_df <- merge(cohort_df, surv_df[, target_cols, drop = FALSE], by = "Group", all.x = TRUE)

          # Fallbacks if columns were missing
          if (!"median" %in% names(cohort_df)) cohort_df$median <- NA
          if (!"0.95LCL" %in% names(cohort_df)) cohort_df$`0.95LCL` <- NA
          if (!"0.95UCL" %in% names(cohort_df)) cohort_df$`0.95UCL` <- NA

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
    }

    # --- 3. COX PROPORTIONAL-HAZARDS ---
    cli::cli_h2("3. Cox Proportional-Hazards")
    formula_str <- if (num_cuts > 0) paste("survival::Surv(time, event) ~ group", cov_part) else paste("survival::Surv(time, event) ~ factor", cov_part)
    model_data <- if (num_cuts > 0) data else object$userdata

    fit_cox <- tryCatch(survival::coxph(stats::as.formula(formula_str), data = model_data), error = function(e) NULL)

    if (is.null(fit_cox)) {
      cli::cli_inform("Could not fit Cox model for best model: convergence failed.")
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

    # --- 4. PROPORTIONAL HAZARDS ASSUMPTION ---
    if (show_ph_test && !is.null(fit_cox)) {
      cli::cli_h2("4. Time-Dependent Diagnostics (Schoenfeld)")

      zph_test <- tryCatch(survival::cox.zph(fit_cox), error = function(e) NULL)

      if (!is.null(zph_test)) {
        zph_table <- zph_test$table
        global_p <- if ("GLOBAL" %in% rownames(zph_table)) zph_table["GLOBAL", "p"] else zph_table[nrow(zph_table), "p"]

        if (!is.na(global_p) && global_p < 0.05) {
          # TIER 2: Hybrid Insight Message
          cli::cli_alert_info("Insight: The hazard ratios appear to fluctuate over time (Global p = {.val {signif(global_p, 3)}}).")
          cli::cli_text("The cut-points successfully separate the subjects, but the relative event risk between these cohorts likely evolves as follow-up time increases.")
        } else if (!is.na(global_p)) {
          # TIER 1: Hybrid Success Message
          cli::cli_alert_success("Passed: The proportional hazards assumption holds across the follow-up period (Global p = {.val {signif(global_p, 3)}}).")
        }
      }
      cat("\n")
    }
  }

  if (plot.it) {
    cli::cli_h2("Model Selection Plot")
    print(plot(object, ...))
    cat("\n")
  }

  # --- 5. PARAMETERS ---
  cli::cli_h2("5. Analysis Parameters")
  params <- object$parameters
  param_bullets <- c(
    "*" = "Search Method: {tools::toTitleCase(params$method %||% \"Unknown\")}",
    "*" = "Predictor: {params$predictor %||% \"Unknown\"}",
    "*" = "Criterion: {params$criterion %||% \"Unknown\"}",
    "*" = "Maximum Cuts Evaluated: {params$max_cuts %||% NA}",
    "*" = "Minimum Group Size (nmin): {params$nmin %||% NA}"
  )
  if (!is.null(params$covariates)) {
    param_bullets <- c(
      param_bullets,
      "*" = "Covariates: {paste(params$covariates, collapse = ', ')}"
    )
  }

  cli::cli_bullets(param_bullets)
  cat("\n")

  invisible(object)
}

#' @rdname find_cutpoint_number
#' @section srrstats compliance:
#' .
#' @srrstats {RE6.0} `plot()` method provided.
#' @srrstats {RE6.2} `plot()` visualises model selection metric vs cuts.
#' @export
plot.find_cutpoint_number_result <- function(x, y, ...) {
  results <- x$results
  criterion_text <- x$parameters$criterion %||% "IC"

  if (is.null(results) || nrow(results) == 0 ||
      all(is.na(results[[criterion_text]]))) {
    cli::cli_inform("Cannot generate plot: no valid IC values found.")
    return(invisible(NULL))
  }

  y_values <- results[[criterion_text]]
  valid_indices <- !is.na(y_values)
  if (sum(valid_indices) == 0) {
    cli::cli_inform("Cannot generate plot: no valid IC values found.")
    return(invisible(NULL))
  }

  plot_data <- results[valid_indices, ]
  y_values <- y_values[valid_indices]

  best_point_idx <- which.min(y_values)
  best_num_cuts <- plot_data$num_cuts[best_point_idx]

  # CRAN SAFEGUARD: Avoid tidyverse global variable warnings by using base mapping
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(
    x = num_cuts,
    y = !!ggplot2::sym(criterion_text)
  )) +
    ggplot2::geom_line(color = "gray50", linewidth = 0.8) +
    ggplot2::geom_point(
      shape = 21,
      size = 3.5,
      fill = "dodgerblue",
      color = "white",
      stroke = 1
    ) +
    ggplot2::geom_point(
      data = function(df) subset(df, num_cuts == best_num_cuts),
      color = "#D55E00", size = 4, shape = 19
    ) +
    ggplot2::scale_x_continuous(breaks = plot_data$num_cuts) +
    ggplot2::labs(
      title = paste("Model Selection by", toupper(criterion_text)),
      subtitle = "Best model (lowest) in orange",
      x = "Number of Cut-points",
      y = toupper(criterion_text)
    ) +
    ggplot2::theme_minimal(base_size = 14) +
    ggplot2::theme(plot.subtitle = ggplot2::element_text(size = 10))

  return(p)
}
