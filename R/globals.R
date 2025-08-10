# This file declares global variables to satisfy R CMD check.
# These are column names created and used within dplyr pipes or ggplot aes().
utils::globalVariables(c(
  # Used in plot.find_cutpoint
  "density",

  # Used in plot.validate_cutpoint_result
  "Value",
  "everything",

  # Used in foreach loop in .systematic_search
  "c1"
))
