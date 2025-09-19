# S3 methods run without errors and produce correct output

    Code
      res_fc
    Message
      
      -- Optimal Cut-point Analysis for Survival Data (Systematic) -------------------
      * Predictor: predictor
      * Criterion: logrank
      * Optimal Log-Rank Statistic: 7.9778
      v Recommended Cut-point(s): 41.51

---

    Code
      summary(res_fc)
    Message
      
      -- Optimal Cut-point Analysis for Survival Data (Systematic) -------------------
      * Predictor: predictor
      * Criterion: logrank
      * Optimal Log-Rank Statistic: 7.9778
      v Recommended Cut-point(s): 41.51
      
      -- Group Counts --
      
    Output
        Group   N Events
      1    G1  74     57
      2    G2 251    176
    Message
      -- Median Survival by Group --
      
    Output
      Call: survfit(formula = survival::Surv(time, event) ~ group, data = data)
      
                 n events median 0.95LCL 0.95UCL
      group=G1  74     57   18.3    13.0    26.8
      group=G2 251    176   24.7    19.6    28.4
    Message
      -- Final Cox Model Summary --
      
    Output
      Call:
      survival::coxph(formula = as.formula(formula_str), data = data)
      
        n= 325, number of events= 233 
      
                 coef exp(coef) se(coef)      z Pr(>|z|)   
      groupG2 -0.4365    0.6463   0.1557 -2.803  0.00507 **
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
              exp(coef) exp(-coef) lower .95 upper .95
      groupG2    0.6463      1.547    0.4763     0.877
      
      Concordance= 0.531  (se = 0.016 )
      Likelihood ratio test= 7.32  on 1 df,   p=0.007
      Wald test            = 7.85  on 1 df,   p=0.005
      Score (logrank) test = 7.98  on 1 df,   p=0.005
      
    Message
      -- Proportional Hazards Assumption Test --
      
    Output
             chisq df    p
      group   1.83  1 0.18
      GLOBAL  1.83  1 0.18
    Message
      -- Analysis Parameters ---------------------------------------------------------
      * Search Method: Systematic
      * Predictor: predictor
      * Number of cuts: 1
      * Minimum group size (nmin): 1

---

    Code
      res_fcn
    Message
      
      -- Optimal Cut-point Number Analysis -------------------------------------------
      Method: systematic
      Criterion: BIC
    Output
       num_cuts     BIC Delta_BIC BIC_Weight                  Evidence  cuts
              0 2219.74      5.27       6.7% Considerably less support    NA
              1 2214.47      0.00      93.3%       Substantial support 41.51
    Message
      v Conclusion: The model with 1 cut-point(s) is the most plausible based on BIC.
      └─ Optimal cuts found at: 41.51
      Hint: Use `summary()` for full model details and `plot()` to visualize this
      table.

---

    Code
      summary(res_fcn)
    Message
      
      -- Optimal Cut-point Number Analysis -------------------------------------------
      Method: systematic
      Criterion: BIC
    Output
       num_cuts     BIC Delta_BIC BIC_Weight                  Evidence  cuts
              0 2219.74      5.27       6.7% Considerably less support    NA
              1 2214.47      0.00      93.3%       Substantial support 41.51
    Message
      v Conclusion: The model with 1 cut-point(s) is the most plausible based on BIC.
      └─ Optimal cuts found at: 41.51
      Hint: Use `summary()` for full model details and `plot()` to visualize this
      table.
      
      -- Details for Best Model ------------------------------------------------------
      The best model found has 1 cut-point(s).
      Cut-point values: 41.51.
      
      -- Group Counts --
      
    Output
        Group   N Events
      1    G1  74     57
      2    G2 251    176
    Message
      -- Median Survival by Group --
      
    Output
      Call: survfit(formula = survival::Surv(time, event) ~ group, data = data)
      
                 n events median 0.95LCL 0.95UCL
      group=G1  74     57   18.3    13.0    26.8
      group=G2 251    176   24.7    19.6    28.4
    Message
      -- Final Cox Proportional-Hazards Model --
      
    Output
      Call:
      survival::coxph(formula = as.formula(formula_str), data = data)
      
        n= 325, number of events= 233 
      
                 coef exp(coef) se(coef)      z Pr(>|z|)   
      groupG2 -0.4365    0.6463   0.1557 -2.803  0.00507 **
      ---
      Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
      
              exp(coef) exp(-coef) lower .95 upper .95
      groupG2    0.6463      1.547    0.4763     0.877
      
      Concordance= 0.531  (se = 0.016 )
      Likelihood ratio test= 7.32  on 1 df,   p=0.007
      Wald test            = 7.85  on 1 df,   p=0.005
      Score (logrank) test = 7.98  on 1 df,   p=0.005
      

