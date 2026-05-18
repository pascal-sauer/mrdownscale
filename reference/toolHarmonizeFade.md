# toolHarmonizeFade

Tool function for creating a harmonized data set with a smooth s-shaped
transition from historic target data to simulated input data.

## Usage

``` r
toolHarmonizeFade(xInput, xTarget, harmonizationPeriod, level = 3)
```

## Arguments

- xInput:

  input data as magpie object

- xTarget:

  target data as magpie object

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- level:

  passed to toolReplaceExpansion

## Value

harmonized data set as magpie object with data from input for years
before the harmonization period, data from target for years after the
harmonization period and a smooth transition in between.

## Author

Jan Philipp Dietrich, Pascal Sauer
