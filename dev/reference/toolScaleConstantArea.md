# toolScaleConstantArea

Scale x to make the sum of all land types is constant over time.

## Usage

``` r
toolScaleConstantArea(x, ..., noteThreshold = 10^-10, warnThreshold = 10^-5)
```

## Arguments

- x:

  a magpie object

- ...:

  not used, will throw an error if supplied

- noteThreshold:

  if the maximum difference between the scaled and the original data is
  greater than this, a note is triggered

- warnThreshold:

  if the maximum difference between the scaled and the original data is
  greater than this, a warning is triggered

## Value

x scaled

## Author

Pascal Sauer
