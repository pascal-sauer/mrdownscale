# toolHandleNegatives

Negative values are set to zero and the remaining values are scaled down
so that the category sums of every cell and year match the target area,
which defaults to the category sums of x before zeroing.

## Usage

``` r
toolHandleNegatives(x, targetArea = NULL, lastScaled = character(0))
```

## Arguments

- x:

  magclass object, may contain negative values

- targetArea:

  optional target area per cell and year, recycled if it has fewer years
  than x, must not be negative

- lastScaled:

  character vector of category names that are only scaled down after all
  other categories have been scaled down to zero, names not present in x
  are ignored

## Value

x without negative values, category sums matching the target area

## Details

Categories listed in lastScaled are only scaled down after all other
categories have been scaled down to zero, so they are protected against
compensation cuts unless the other categories cannot absorb the excess
area alone. With the default lastScaled = character(0) all categories
are scaled down proportionally together.

Intended to be applied to a group of land cover categories (e.g. forest
plus other land), so that negative values in one category are
compensated by the others without changing the total group area.

## Author

Pascal Sauer
