# calcLandTransitions

Ex-Post estimation of land use transitions based on land use state
information

## Usage

``` r
calcLandTransitions(
  outputFormat,
  input,
  harmonizationPeriod,
  yearsSubset,
  harmonization,
  downscaling,
  gross
)
```

## Arguments

- outputFormat:

  format in which the outputs should be prepared. Currently, only "ESM"
  for earth system model compatible input data is available.

- input:

  name of an input dataset, currently only "magpie"

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- yearsSubset:

  vector of years to keep in the output dataset

- harmonization:

  name of harmonization method, see
  [`toolGetHarmonizer`](toolGetHarmonizer.md)

- downscaling:

  name of downscaling method, currently only "magpieClassic"

- gross:

  either boolean or a magpie object containing bidirectional transition
  shares relative to the area of the involved land pools (transition
  divided by the area of the land pool in the "from" sub dimension). If
  set to FALSE only net transitions will be returned. If set to TRUE an
  internal gross transition estimate based on average gross transitions
  in LUH2 in the period from 1995 to 2015 will be used.

## Value

land use transition data

## Author

Jan Philipp Dietrich, Pascal Sauer
