# calcNonlandHighRes

Calculate a high resolution dataset from the low resolution input
dataset and high resolution data.

## Usage

``` r
calcNonlandHighRes(
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

  name of an input dataset, currently only "magpie"

- target:

  name of a target dataset, see [`calcLandTarget`](calcLandTarget.md)
  for available target datasets

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

downscaled nonland data

## Details

Wood harvest area is disaggregated using the maximum possible harvest
per year which is based on the downscaled land data. Bioh is
disaggregated using the just disaggregated wood harvest area as weight.
Fertilizer is in kg ha-1 yr-1, so we simply use the low res/region value
for each cell corresponding to that region. Harvest weight type is
disaggregated using the nonland target data in the first year of the
harmonization period as weight.

## Author

Pascal Sauer

## Examples

``` r
if (FALSE) { # \dontrun{
  calcOutput("NonlandHighRes", input = "magpie", target = "luh2mod",
             harmonizationPeriod = c(2015, 2050), yearsSubset = 2015:2100,
             harmonization = "fade", downscaling = "magpieClassic")
} # }
```
