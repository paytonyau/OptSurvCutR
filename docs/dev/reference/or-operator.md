# Infix operator for NULL default

Infix operator for NULL default

## Usage

``` r
a %||% b
```

## Arguments

- a:

  The value to check.

- b:

  The default value if \`a\` is \`NULL\`.

## Value

\`a\` if not \`NULL\`, otherwise \`b\`.

## srrstats compliance

.

## Examples

``` r
x <- NULL
y <- 5
x %||% y # Returns 5
#> [1] 5

z <- 10
z %||% y # Returns 10
#> [1] 10
```
