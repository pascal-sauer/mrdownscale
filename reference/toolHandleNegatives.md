# toolHandleNegatives

Negative values are set to zero and the remaining values are scaled down
so that the category sums of every cell and year match the target area,
which defaults to the category sums of x before zeroing.

## Usage

``` r
toolHandleNegatives(x, targetArea = NULL)
```

## Arguments

- x:

  magclass object, may contain negative values

- targetArea:

  optional target area per cell and year, recycled if it has fewer years
  than x, must not be negative

## Value

x without negative values, category sums matching the target area

## Details

Intended to be applied to a group of land cover categories (e.g. forest
plus other land), so that negative values in one category are
compensated by the others without changing the total group area.

## Author

Pascal Sauer
