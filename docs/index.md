# OptSurvCutR: Validated Cut-point Selection for Survival Analysis

![OptSurvCutR logo](reference/figures/logo.png)

`OptSurvCutR` (**Opt**imal **Surv**ival **Cut**-points **R**) provides a
rigorous, reproducible **three-step workflow** for discovering the
optimal number and location of cut-points in time-to-event (survival)
data. Designed for continuous predictors (e.g., gene expression, virome
abundance, biomarkers), it moves beyond arbitrary median splits to fully
**data-driven stratification**.

![OptSurvCutR Graphical
Abstract](reference/figures/graphical_abstract.jpg)

OptSurvCutR Graphical Abstract

## Why OptSurvCutR?

| Feature | Benefit |
|----|----|
| Optimal number of cuts | Uses AIC, AICc, or BIC to select 0–k cut-points |
| Flexible search | Systematic grid or genetic algorithm (`rgenoud`) |
| Covariate adjustment | Control for confounders during cut-point discovery |
| Bootstrap validation | 95% confidence intervals for cut-point stability |
| Publication-ready plots | Kaplan–Meier, optimisation curves, Schoenfeld residual diagnostics |

## Installation

You can install the development version of `OptSurvCutR` from GitHub.
Note that the genetic algorithm (`method = "genetic"`) requires the
`rgenoud` package, which should be installed separately from CRAN if you
plan to use it.

``` r
# Core dependencies
install.packages(c("remotes", "survival"))

# Optional but highly recommended for >2 cuts
install.packages("rgenoud")

# Install the development version from GitHub
remotes::install_github("paytonyau/OptSurvCutR")
```

## Example: Quick Workflow with CRC Virome Data

``` r
library(OptSurvCutR)
library(survival)

data("crc_virome")

crc <- crc_virome %>%
  select(time = time_months, status, Enterovirus) %>%
  na.omit()

# Step 1: Determine optimal number of cut-points
num_cuts <- find_cutpoint_number(
  data = crc, predictor = "Enterovirus",
  outcome_time = "time", outcome_event = "status",
  max_cuts = 2, nmin = 0.15, seed = 42
)
print(num_cuts)   # BIC suggests 2 cut-points

# Step 2: Find the precise cut-point locations
cuts <- find_cutpoint(
  data = crc, predictor = "Enterovirus",
  outcome_time = "time", outcome_event = "status",
  num_cuts = 2, method = "systematic", nmin = 0.15
)

# Step 3: Validate stability with bootstrap
val <- validate_cutpoint(cuts, num_replicates = 200, seed = 456)
summary(val)

# Step 4: Visualise results
plot(cuts, type = "outcome")        # Kaplan–Meier curves
plot(cuts, type = "distribution")  # Predictor + cuts
plot_schoenfeld(cuts)               # Proportional hazards check
plot(val)                           # Bootstrap density + 95% CI
```

## Workflow Summary

OptSurvCutR provides a three-step workflow for cut-point analysis:

1.  [`find_cutpoint_number()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint_number.md):
    – selects the statistically optimal number of cut-points using
    information criteria
2.  [`find_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/find_cutpoint.md):
    locates exact cut-point values (systematic or genetic search)
3.  [`validate_cutpoint()`](https://paytonyau.github.io/OptSurvCutR/reference/validate_cutpoint.md):
    assesses stability via bootstrapping with 95% confidence intervals
    \## Resources

- **Vignettes**: `browseVignettes("OptSurvCutR")`
- **Package Website**: <https://paytonyau.github.io/OptSurvCutR/>
- **Manuscript**: Yau, Payton T. O. “OptSurvCutR: Validated Cut-point
  Selection for Survival Analysis.” bioRxiv preprint, posted October
  10, 2025. <https://doi.org/10.1101/2025.10.08.681246>.

## Citation

``` r
@article{yau2025optsurvcutr,
  author  = {Yau, Payton T. O.},
  title   = {OptSurvCutR: Validated Cut-point Selection for Survival Analysis},
  year    = {2025},
  doi     = {10.1101/2025.10.08.681246},
  publisher = {Cold Spring Harbor Laboratory},
  journal = {bioRxiv},
  url     = {https://www.biorxiv.org/content/10.1101/2025.10.08.681246}
}
```

A JOSS submission is planned post-rOpenSci review.

## Support OptSurvCutR

If `OptSurvCutR` helps your research, consider buying me a coffee — it
directly supports ongoing maintenance with no dedicated funding.

[![Buy Me A
Coffee](https://img.buymeacoffee.com/button-api/?text=Buy%20me%20a%20coffee&emoji=&slug=payton.yau&button_colour=FFDD00&font_colour=000000&font_family=Poppins&outline_colour=000000&coffee_colour=ffffff)](https://buymeacoffee.com/payton.yau)

## License

Licensed under the GPL-3 License.

## Contact

Questions, suggestions, or issues? Please open a ticket:
<https://github.com/paytonyau/OptSurvCutR/issues>
