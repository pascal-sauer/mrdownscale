# toolCheckWoodHarvestArea

Check wood harvest area is not exceeding land area of the corresponding
type divided by timestep length. Also, check that primf and primn are
reduced by at least as much as they were harvested.

## Usage

``` r
toolCheckWoodHarvestArea(harvest, land, endOfHistory)
```

## Arguments

- harvest:

  magpie object with exactly the following categories: paste0(c("primf",
  "secyf", "secmf", "pltns", "primn", "secnf"), "\_wood_harvest_area")

- land:

  magpie object with at least the following categories: c("primf",
  "secdf", "pltns", "primn", "secdn")

- endOfHistory:

  The last year considered part of the historical period, will check and
  report consistency separately for history and after

## Author

Pascal Sauer
