# This file declares global variables to satisfy R CMD check.
# These are column names created and used within dplyr pipes or ggplot aes().
utils::globalVariables(c(
  # From previous versions
  "density", "Value", "everything", "c1",

  # New variables from plot_effect_size and plot_waterfall
  "Cut1", "HR", "HR_low", "HR_up", "OR", "OR_low", "OR_up",
  "classified_group", "is_correct", "outcome", "patient_id"
))
