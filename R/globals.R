# This file declares global variables to satisfy R CMD check.
# These are column names created and used within dplyr pipes or ggplot aes().
utils::globalVariables(c(
  # From previous versions
  "density", "Value", "everything", "c1",

  # New variables from plot_effect_size and plot_waterfall
  "Cut1", "HR", "HR_low", "HR_up", "OR", "OR_low", "OR_up",
  "classified_group", "is_correct", "outcome", "patient_id",

  # for plotting_functions.R
  "cut1",

  # --- ADDED FROM DEBUG ---
  "stat",               # Used in plot_optimization_curve
  "factor",             # Used in plot.find_cutpoint
  "original_cut",       # Used in plot.validate_cutpoint
  "ci_lower",           # Used in plot.validate_cutpoint
  "ci_upper"            # Used in plot.validate_cutpoint
))
