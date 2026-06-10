# Diagnostic Plot of Schoenfeld Residuals

High-tier multi-panel diagnostic dashboard tracking the proportional
hazards assumption with custom facets per risk cohort stratum.

## Usage

``` r
plot_cutpoint_residuals(x, ...)
```

## Arguments

- x:

  A `find_cutpoint` result object.

- ...:

  Unused optional arguments.

## Value

A publication-ready `ggplot` canvas frame, or `NULL` if the fit fails.

## Examples

``` r
mock_df <- data.frame(time = 1:30, event = rep(c(0, 1), 15), factor = rnorm(30))
res <- find_cutpoint(
  mock_df, "factor", "time", "event",
  num_cuts = 1, method = "systematic", quiet = TRUE
)
p <- plot_cutpoint_residuals(res)
#> No valid cut-points mapped; diagnostics skipped.
```
