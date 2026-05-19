# calcLandTargetLowRes

Aggregate target land data to the spatial resolution of the input data
in preparation for harmonization.

## Usage

``` r
calcLandTargetLowRes(input, target, endOfHistory)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  name of a target dataset, see [`calcLandTarget`](calcLandTarget.md)
  for available target datasets

- endOfHistory:

  years later than this are not returned

## Value

low resolution target land data

## Author

Pascal Sauer
