# calcStatesNC

Prepare data to be written as LUH-style states.nc file. Call this via
calcOutput in a full function, and set calcOutput's file argument to a
.nc file path.

## Usage

``` r
calcStatesNC(
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

data prepared to be written as a LUH-style states.nc file

## Author

Pascal Sauer, Jan Philipp Dietrich
