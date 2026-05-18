# toolAggregateCropland

Aggregate variables like c3ann_irrigated and c3ann_rainfed to just
c3ann.

## Usage

``` r
toolAggregateCropland(
  land,
  cropTypes = c("c3ann", "c4ann", "c3per", "c4per", "c3nfx"),
  ...,
  keepOthers = TRUE
)
```

## Arguments

- land:

  magpie object with variables starting with cropTypes in dim 3

- cropTypes:

  character vector, the variables to aggregate to

- ...:

  reserved for future expansion

- keepOthers:

  logical, if FALSE cropTypes will be the only variables in output

## Value

land with crop data aggregated to cropTypes

## Author

Pascal Sauer
