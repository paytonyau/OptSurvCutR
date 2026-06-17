# Plot Optimisation Curve or Surface from Search

Plots the metric landscape evaluated across coordinates. Maps a 1D
optimisation line for 1-cut systematic setups, or a 2D topographic grid
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
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  # Build clean simulation objects using a 2-cut systematic search
  # to populate the 2D grid matrix log array
  set.seed(123)
  mock_df <- data.frame(
    time   = c(runif(10, 50, 100), runif(10, 20, 60), runif(10, 5, 25)),
    event  = rep(1, 30),
    factor = c(rnorm(10, 2, 0.2), rnorm(10, 7, 0.2), rnorm(10, 15, 0.2))
  )
  res <- find_cutpoint(
    mock_df, "factor", "time", "event",
    num_cuts = 2, method = "systematic", quiet = TRUE, nmin = 3
  )
  p <- plot_optimisation_curve(res)
}
```
