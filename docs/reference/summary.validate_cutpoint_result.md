# Summary of Bootstrap Stability Validation

Calculates precision and boundary validity widths to route the
validation runs into an automated multi-tier clinical suitability score.

## Usage

``` r
# S3 method for class 'validate_cutpoint_result'
summary(
  object,
  show_descriptives = TRUE,
  show_ci = TRUE,
  show_params = TRUE,
  plot.it = FALSE,
  ...
)
```

## Arguments

- object:

  A `validate_cutpoint_result` object.

- show_descriptives:

  Logical. Show complete distribution descriptives?

- show_ci:

  Logical. Print 95% Confidence Interval boundaries?

- show_params:

  Logical. Display execution tracking parameters?

- plot.it:

  Logical. If `TRUE`, automatically prints the sampling line chart.

- ...:

  Additional arguments passed down to downstream rendering pipelines.

## Value

The validation object invisibly.

## srrstats compliance

.
