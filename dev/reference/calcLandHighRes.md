# calcLandHighRes

This function performs the downscaling: It calculates a high resolution
dataset from the low resolution input dataset and the high resolution
target dataset using the given downscaling method.

## Usage

``` r
calcLandHighRes(
  input,
  target,
  harmonizationPeriod,
  yearsSubset,
  harmonization,
  downscaling
)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  name of a target dataset

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

downscaled land use data

## Author

Jan Philipp Dietrich, Pascal Sauer
