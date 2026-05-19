# calcNonlandHarmonized

Harmonize nonland input data to target data using the specified method,
checking data for consistency before returning.

## Usage

``` r
calcNonlandHarmonized(input, target, harmonizationPeriod, harmonization)
```

## Arguments

- input:

  name of an input dataset, currently only "magpie"

- target:

  name of the target dataset

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- harmonization:

  harmonization method, see [`toolGetHarmonizer`](toolGetHarmonizer.md)
  for available methods

## Value

harmonized nonland data

## Details

Wood harvest biomass (bioh) is adapted to the harmonized wood harvest
area by calculating kg C per mega hectare for input and target data and
harmonizing it. This is then multiplied by the harmonized wood harvest
area and scaled so the total harmonized bioh is reached. Harmonize
absolute fertilizer in Tg yr-1, then convert to fertilizer rate in kg
ha-1 yr-1.

## Author

Pascal Sauer
