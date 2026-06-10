# Custom Clinical Theme for OptSurvCutR

A clean, publication-ready `ggplot2` theme used as the standard
aesthetic for all OptSurvCutR plots. Features subtle light-grey
background gridlines for high-precision tracing.

## Usage

``` r
theme_optsurv(base_size = 14)
```

## Arguments

- base_size:

  Base font size, default is 14.

## Value

A `ggplot2` theme object containing customized layout parameters.

## Examples

``` r
library(ggplot2)
mock_df <- data.frame(x = 1:10, y = 1:10)
ggplot(mock_df, aes(x, y)) +
  geom_point() +
  theme_optsurv()
```
