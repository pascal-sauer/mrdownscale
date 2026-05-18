# calcNonlandInput

Prepare the nonland input data for category mapping, checking data for
consistency before returning.

## Usage

``` r
calcNonlandInput(input)
```

## Arguments

- input:

  name of an input dataset, currently only "magpie"

## Value

nonland input data

## Details

All "Land" functions deal with area data, as opposed to "Nonland"
functions which deal with non-area data such as the amount of applied
fertilizer. These are treated differently, because for area data other
constraints apply, e.g. the total area must be constant over time.
Fertilizer on regional level is disaggregated to cluster level using
cropland as weight.

## Author

Pascal Sauer
