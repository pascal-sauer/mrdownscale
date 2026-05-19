# calcNonlandTarget

Prepare the high resolution nonland target dataset for harmonization and
downscaling, checking data for consistency before returning.

## Usage

``` r
calcNonlandTarget(target, endOfHistory)
```

## Arguments

- target:

  name of a target dataset, see [`calcLandTarget`](calcLandTarget.md)
  for available target datasets

- endOfHistory:

  years later than this are not returned

## Value

nonland target data

## Author

Pascal Sauer
