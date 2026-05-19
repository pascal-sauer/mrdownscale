# fullDOWNSCALEDMAGPIE

Run the pipeline to generate harmonized and downscaled MAgPIE data using
landuseinit as a reference dataset. Write output in the format of
avl_land_t_0.5.mz, full report on consistency checks is printed and
written to report.log.

## Usage

``` r
fullDOWNSCALEDMAGPIE(
  rev = numeric_version("0"),
  harmonizationPeriod = c(2015, 2050),
  target = "landuseinit",
  downscaling = "magpieClassic",
  harmonization = "fade"
)
```

## Arguments

- rev:

  revision number of the data. If not provided the current date will be
  used instead. When called via madrat::retrieveData rev will be
  converted to numeric_version.

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- target:

  Name of dataset to be used as harmonization target and downscaling
  reference

- downscaling:

  name of downscaling method, currently only "magpieClassic"

- harmonization:

  name of harmonization method, see
  [`toolGetHarmonizer`](toolGetHarmonizer.md)

## Author

Pascal Sauer
