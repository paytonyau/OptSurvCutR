# Coverage: summary.find_cutpoint_number_result – no valid IC

    Code
      summary(obj)
    Message
      
      -- Optimal Cut-point Number Analysis (Systematic) ------------------------------
      Cannot summarise: no valid model was found.

# Coverage: plot.find_cutpoint_number_result – all IC NA

    Code
      plot(obj)
    Message
      Cannot generate plot: no valid IC values found.

# summary handles pathological IC (mix of NA/Inf)

    Code
      summary(obj)
    Message
      
      -- Optimal Cut-point Number Analysis (Systematic) ------------------------------
      Cannot summarise: no valid model was found.

# S3 methods handle missing parameters gracefully

    Code
      print(obj)
    Message
      
      -- Optimal Cut-point Number Analysis -------------------------------------------
      Method: Unknown
      Criterion: IC
      No optimal model could be determined.

---

    Code
      summary(obj)
    Message
      
      -- Optimal Cut-point Number Analysis (Unknown) ---------------------------------
      Cannot summarise: no valid model was found.

---

    Code
      plot(obj)
    Message
      Cannot generate plot: no valid IC values found.

