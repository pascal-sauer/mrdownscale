# calcNonlandTargetLowRes

Aggregate target nonland data to the spatial resolution of the input
data in preparation for harmonization. Fertilizer is converted to Tg
yr-1, then aggregated and converted back to kg ha-1 yr-1.

## Usage

``` r
calcNonlandTargetLowRes(input, target, endOfHistory)
```

## Arguments

- input:

  name of an input dataset, currently only "magpie"

- target:

  name of a target dataset

- endOfHistory:

  years later than this are not returned

## Value

low resolution target nonland data

## Author

Pascal Sauer
