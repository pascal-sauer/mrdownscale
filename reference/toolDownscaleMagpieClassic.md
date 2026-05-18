# toolDownscaleMagpieClassic

classic MAgPIE downscaling method using
[`luscale::interpolate2`](https://rdrr.io/pkg/luscale/man/interpolate2.html)
as downscaling function.

## Usage

``` r
toolDownscaleMagpieClassic(x, xTarget, xTargetLowRes, mapping)
```

## Arguments

- x:

  magclass containing land to be downscaled

- xTarget:

  magclass target land use dataset for initialization year

- xTargetLowRes:

  magclass target land use dataset for initialization year in low
  resolution (like x)

- mapping:

  mapping between `x` and `xTarget`

## Value

downscaled land use dataset

## Author

Jan Philipp Dietrich, Pascal Sauer
