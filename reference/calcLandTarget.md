# calcLandTarget

Prepare the high resolution target land use dataset for harmonization
and downscaling, checking data for consistency before returning.

## Usage

``` r
calcLandTarget(target, endOfHistory)
```

## Arguments

- target:

  name of the target dataset, one of luh2, luh2mod, luh3 luh2mod/luh3
  will split secdf into pltns and secdf

- endOfHistory:

  years later than this are not returned

## Value

land target data

## Author

Pascal Sauer
