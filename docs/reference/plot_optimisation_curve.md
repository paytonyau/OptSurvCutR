# Plot Optimisation Curve or Surface from Search

Plots the metric landscape evaluated across coordinates. Maps a 1D
optimization line for 1-cut systematic setups, or a 2D topographic grid
profile for 2-cut layouts.

## Usage

``` r
plot_optimisation_curve(cutpoint_result, ...)
```

## Arguments

- cutpoint_result:

  A `find_cutpoint` object.

- ...:

  Unused dots.

## Value

A valid `ggplot` object detailing evaluation statistics vs coordinates.

## srrstats compliance

.

## Examples

``` r
mock_df <- data.frame(time = 1:20, event = rep(c(0, 1), 10), factor = rnorm(20, 10, 2))
res <- find_cutpoint(mock_df, "factor", "time", "event", num_cuts = 1, method = "systematic", quiet = TRUE)
p <- plot_optimisation_curve(res)
#> Error in plot_optimisation_curve(res): The results object must contain a valid grid log array in `all_stats`.
```
