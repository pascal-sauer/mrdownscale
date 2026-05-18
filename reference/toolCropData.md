# toolCropData

For each crop type (c3ann, c3nfx, c3per, c4ann, c4per) calculate
irrigation share, 1st and 2nd generation biofuel shares and share of
cell area (corresponding to LUH3 variables c3ann, irrig_c3ann,
cpbf1_c3ann, cpbf2_c3ann and analogously for the other crop types).

## Usage

``` r
toolCropData(landHighRes, cellArea)
```

## Arguments

- landHighRes:

  high resolution land use data as magclass object

- cellArea:

  corresponding magclass object containing cell area in Mha

## Value

crop data as magclass object with variables irrig\_\*, cpbf1\_\*,
cpbf2\_\*, \* for all 5 crop types

## Author

Pascal Sauer
