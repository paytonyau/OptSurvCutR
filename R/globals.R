# This file declares global variables to satisfy R CMD check.
# These are column names created and used within dplyr pipes or ggplot aes().
utils::globalVariables(c(
  ".", "Cut1", "HR", "HR_low", "HR_up", "OR", "OR_low", "OR_up",
  "is_correct", "patient_id", "classified_group", "outcome", "factor",
  "abundance", "abundance_group", "time", "status", "sample_id", "status_text",
  "NAME", "predictor", "Group", "Group_Status", "p_value", "sd_AIC",
  "mean_AIC", "Delta_AIC", "Akaike_Weight", "Delta_BIC", "BIC_Weight",
  "Value", "Cutpoint", "specificity", "sensitivity"
))
