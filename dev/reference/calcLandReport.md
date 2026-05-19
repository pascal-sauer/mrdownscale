# calcLandReport

Convert the downscaled land use data to the format required by the given
project.

## Usage

``` r
calcLandReport(
  outputFormat,
  input,
  harmonizationPeriod,
  yearsSubset,
  harmonization,
  downscaling
)
```

## Arguments

- outputFormat:

  format in which the outputs should be prepared. Options: ESM,
  ScenarioMIP, downscaledmagpie

- input:

  name of an input dataset, options: "magpie", "witch"

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

## Value

land use data

## Author

Pascal Sauer
