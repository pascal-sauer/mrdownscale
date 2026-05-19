# fullESM

Run the pipeline to generate harmonized and downscaled data to report
for the RESCUE, OptimESM and other projects where ESM compatible land
use inputs are required. Write .nc files, print full report on
consistency checks and write it to report.log.

## Usage

``` r
fullESM(
  rev = numeric_version("0"),
  input = "magpie",
  scenario = "",
  harmonizationPeriod = c(2015, 2050),
  yearsSubset = 2015:2100,
  harmonization = "fade",
  downscaling = "magpieClassic",
  compression = 2,
  progress = TRUE
)
```

## Arguments

- rev:

  revision number of the data. If not provided the current date will be
  used instead. When called via madrat::retrieveData rev will be
  converted to numeric_version.

- input:

  name of an input dataset, options: "magpie", "witch"

- scenario:

  scenario name to be included in filenames

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

- compression:

  compression level of the resulting .nc files, possible values are
  integers from 1-9, 1 = fastest, 9 = best compression

- progress:

  boolean defining whether progress should be printed

## Author

Pascal Sauer, Jan Philipp Dietrich
