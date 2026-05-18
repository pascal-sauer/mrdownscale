# toolHarmonizeFadeForest

Harmonize using [`toolHarmonizeFade`](toolHarmonizeFade.md), then sum up
primf + secdf, then disaggregate again in a way that minimizes primf to
secdf conversion, without exceeding linear extrapolation of primf.

## Usage

``` r
toolHarmonizeFadeForest(xInput, xTarget, harmonizationPeriod)
```

## Arguments

- xInput:

  input data as magpie object

- xTarget:

  target data as magpie object

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

## Value

harmonized data set as magpie object with data from input for years
before the harmonization period, data from target for years after the
harmonization period and a smooth transition in between.

## Author

Pascal Sauer
