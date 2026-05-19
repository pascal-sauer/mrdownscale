# toolMaxHarvestPerYear

Calculate the maximum possible harvest area per year based on the given
land data.

## Usage

``` r
toolMaxHarvestPerYear(land, split = TRUE, timestepAdjust = TRUE)
```

## Arguments

- land:

  magpie object with at least the following categories: c("primf",
  "secdf", "pltns", "primn", "secdn")

- split:

  if TRUE: split secdf to secyf and secdmf, rename secdf to secnf, and
  add dim "wood_harvest_area"

- timestepAdjust:

  if TRUE: divide values for primary land by timestep length. This makes
  sense, because once primary land has been harvested, it is converted
  to secondary land and thus cannot be harvested again. Might introduce
  unintended spikes when timestep length changes.

## Value

magpie object with the maximum possible yearly wood harvest area

## Author

Pascal Sauer
