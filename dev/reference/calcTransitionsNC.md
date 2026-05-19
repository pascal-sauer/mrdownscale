# calcTransitionsNC

Prepared data to be written as a LUH-style transitions.nc file

## Usage

``` r
calcTransitionsNC(
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

  options: ESM, ScenarioMIP

- input:

  name of an input dataset, options: "magpie", "witch"

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- yearsSubset:

  remove years from the returned data which are not in yearsSubset

- harmonization:

  name of harmonization method, see
  [`toolGetHarmonizer`](toolGetHarmonizer.md)

- downscaling:

  name of downscaling method, currently only "magpieClassic"

## Value

data prepared to be written as a LUH-style transitions.nc file

## Author

Pascal Sauer, Jan Philipp Dietrich
