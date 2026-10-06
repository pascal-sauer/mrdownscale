# toolReportAreaDeviation

Reports, per category group and per variable, how much the corrected
output deviates from the raw input data: the positive and the negative
deviations (xOut - xRaw) accumulated separately over all cells and
timesteps, as mean per timestep in Mha and as a percentage of the output
area accumulated over the same timesteps, plus the net change as the sum
of both columns, in six columns: =

## Usage

``` r
toolReportAreaDeviation(xRaw, xOut, groups)
```

## Arguments

- xRaw:

  harmonized data before corrections, as magpie object, with the years
  of xOut after the harmonization year

- xOut:

  harmonized data after corrections, as magpie object, including the
  harmonization year and the years before it

- groups:

  Named list of character vectors with the category groups to be
  reported separately, named as they should appear in the report; every
  category must be in one of the groups or be "urban"

## Author

Pascal Sauer
