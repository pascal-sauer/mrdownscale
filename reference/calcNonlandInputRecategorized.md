# calcNonlandInputRecategorized

Harmonize categories by mapping nonland input data categories to the
categories of the nonland target dataset. See
[`calcLandInputRecategorized`](calcLandInputRecategorized.md) for an
explanation of the mapping procedure.

## Usage

``` r
calcNonlandInputRecategorized(
  input,
  target,
  youngShareWoodHarvestArea = 0.95,
  youngShareWoodHarvestWeight = 0.5
)
```

## Arguments

- input:

  name of an input dataset, currently only "magpie"

- target:

  name of a target dataset, see [`calcLandTarget`](calcLandTarget.md)
  for available target datasets

- youngShareWoodHarvestArea:

  share of wood harvest area taken from young (instead of mature)
  secondary forest; default value is based on LUH value from 2014; used
  to disaggregate wood harvest area from secondary forest to secondary
  young and mature forest

- youngShareWoodHarvestWeight:

  analogue to youngShareWoodHarvestArea for wood harvest weight instead
  of area

## Value

nonland data with target categories

## Details

Report and discard wood harvest area if there is zero wood harvest
(bioh) or vice versa.

## Author

Pascal Sauer
