# calcNonlandReport

Convert the downscaled nonland data to the format required by the given
project.

## Usage

``` r
calcNonlandReport(
  outputFormat,
  harmonizationPeriod,
  yearsSubset,
  harmonization,
  downscaling
)
```

## Arguments

- outputFormat:

  options: ESM, ScenarioMIP

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

nonland data

## Author

Pascal Sauer

## Examples

``` r
if (FALSE) { # \dontrun{
  calcOutput("NonlandReport", outputFormat = "ESM",
             harmonizationPeriod = c(2015, 2050), yearsSubset = 2015:2100,
             harmonization = "fade", downscaling = "magpieClassic")
} # }
```
