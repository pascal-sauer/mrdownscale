# calcLandHarmonized

This function computes a version of the chosen land input data in the
resolution of the land input data but harmonized to the land use
information of the chosen land target data set (harmonized categories as
well as harmonized transition from historic target data to simulated
input data).

## Usage

``` r
calcLandHarmonized(input, target, harmonizationPeriod, harmonization)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  name of the land target source to be used

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- harmonization:

  name of harmonization method, see
  [`toolGetHarmonizer`](toolGetHarmonizer.md)

## Author

Pascal Sauer, Jan Philipp Dietrich
