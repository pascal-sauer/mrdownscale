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

  Two years A and B, integer values. Up until (including) A the target
  data is returned unchanged. For harmonization != "absoluteChanges" and
  years between A and B (exclusive, A \< x \< B) harmonized data is
  returned. For years \>= B return input data. For harmonization =
  "absoluteChanges", A must equal B. For years \> A return harmonized
  data (constant offset to input data).

- harmonization:

  harmonization method, see [`toolGetHarmonizer`](toolGetHarmonizer.md)
  for available methods

## Value

harmonized nonland data

## Details

Wood harvest biomass (bioh) is adapted to the harmonized wood harvest
area by calculating kg C per mega hectare for input and target data and
harmonizing it. This is then multiplied by the harmonized wood harvest
area and scaled so the total harmonized bioh is reached. bioh of
clusters whose harmonized harvest area is zero is redistributed to
clusters with harvest area in the same region, falling back to the
global pool (toolRedistributeBioh). Fertilizer is harmonized in Tg yr-1,
and converted to fertilizer rate (kg ha-1 yr-1) at the very end. For
harmonization = "absoluteChanges" harmonize at regional resolution by
starting at the anchor (target data in harmonization year) and applying
input deltas (input - input in harmonization year). Regional data is
then distributed to clusters (weight: harmonized cropland, see
toolProjectFertilizerTg). The rates are then rescaled per region, crop
category and year so that each region's fertilizer total from the
harmonization year on keeps a constant offset to input data (parallel
year-to-year changes), redistributing amounts lost where zero clamping
or missing cropland removed fertilizer to the other cells of the same
region, falling back to the global factor exactly in those crop group
and year tables where a region lost all its cropland and can no longer
hold any fertilizer, and finally capped at the plausibility threshold
(toolCapFertilizerKgPerHa).

## Author

Pascal Sauer
