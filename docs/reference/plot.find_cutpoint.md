# Master S3 Plot Router for find_cutpoint

Unified plotting dispatch network routing to clinical survival curves,
predictor density distributions, hazard ratio forest charts,
multi-dimensional objective surfaces, or conditional landmark
stratification assets.

## Usage

``` r
# S3 method for class 'find_cutpoint'
plot(
  x,
  type = c("outcome", "distribution", "forest", "surface", "trajectory", "diagnostic",
    "landmark"),
  return_data = FALSE,
  landmark = NULL,
  ...
)
```

## Arguments

- x:

  A `find_cutpoint` result object.

- type:

  Plot framework type: `"outcome"`, `"distribution"`, `"forest"`,
  `"surface"`, `"trajectory"`, `"diagnostic"`, or `"landmark"`.

- return_data:

  Logical. If `TRUE`, exits the router early and returns the assigned
  underlying data frame template.

- landmark:

  Numeric. The operational milestone timestamp used if
  `type = "landmark"`.

- ...:

  Additional arguments passed down to downstream rendering pipelines.

## Value

A `ggplot` canvas object, a multi-panel `patchwork` collection, or a
`data.frame` if `return_data = TRUE`.

## srrstats compliance

.

## Examples

``` r
# Build clean local simulation objects with an explicit survival risk split
set.seed(123)
mock_df <- data.frame(
  time   = c(runif(15, 50, 100), runif(15, 5, 25)),
  event  = rep(1, 30),
  factor = c(rnorm(15, 5, 0.5), rnorm(15, 15, 0.5))
)
res <- find_cutpoint(
  mock_df, "factor", "time", "event",
  num_cuts = 1, method = "systematic", quiet = TRUE, nmin = 3
)
p <- plot(res, type = "distribution")
```
