# toolFillYearsSpatRaster

Fill data for missing years using linear interpolation.

## Usage

``` r
toolFillYearsSpatRaster(x, years = NULL)
```

## Arguments

- x:

  SpatRaster with years in layer names

- years:

  data for these years will be added if they are not already present, if
  years is NULL all years between the first and last year in x will be
  filled

## Value

SpatRaster with filled years

## Author

Pascal Sauer
