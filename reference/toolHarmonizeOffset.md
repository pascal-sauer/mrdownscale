# toolHarmonizeOffset

Harmonize datasets using the offset method as implemented in
mip::harmonize, originally implemented in Python as part of aneris.

## Usage

``` r
toolHarmonizeOffset(xInput, xTarget, harmonizationPeriod)
```

## Arguments

- xInput:

  magpie object to fade to, usually model projections

- xTarget:

  magpie object to fade from, usually historical data, this argument is
  called "target" for consistency with argument names of other
  functions, it is somewhat misleading here, because data after the
  harmonizationPeriod will not be taken from target, but from input

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

## Value

magpie object with harmonized data

## Author

Pascal Sauer
