# Plot Landmark Stratification Curves

Computes and visualises conditional survival probabilities for patients
who survive up to a specified landmark milestone.

## Usage

``` r
plot_landmark_stratification(x, landmark = NULL, ...)
```

## Arguments

- x:

  A `find_cutpoint` result object.

- landmark:

  Numeric value indicating the landmark milestone. If NULL, defaults to
  20% of maximum follow-up data timelines.

- ...:

  Additional arguments passed down to downstream rendering pipelines.

## Value

A `ggplot` template or `survminer` survival curve asset mapping layout.

## Examples

``` r
if (requireNamespace("survival", quietly = TRUE)) {
  library(survival)
  mock_df <- data.frame(
    time   = runif(25, 5, 50),
    event  = sample(c(0, 1), 25, replace = TRUE),
    factor = rnorm(25, 10, 2)
  )
  res <- find_cutpoint(
    mock_df, "factor", "time", "event",
    num_cuts = 1, method = "systematic", quiet = TRUE, nmin = 3
  )
  p <- plot_landmark_stratification(res, landmark = 10)
}
#> ℹ Generating Landmark Survival Curve for survivors remaining at time milestone: 10
```
