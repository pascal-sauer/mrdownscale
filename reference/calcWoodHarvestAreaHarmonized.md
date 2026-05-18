# calcWoodHarvestAreaHarmonized

Harmonize wood harvest area based on harmonized land data. First, wood
harvest area is harmonized just like land data. Then, to ensure
consistency, primary harvest is converted to seondary harvest (or vice
versa) to match the primary land reduction from the harmonized land
data. If secondary harvest exceeds the available secondary land, the
excess is shifted to the other secondary land (secdf to secdn and vice
versa). Remaining excess harvest area is reported.

## Usage

``` r
calcWoodHarvestAreaHarmonized(
  input,
  target,
  harmonizationPeriod,
  harmonization
)
```

## Arguments

- input:

  name of an input dataset, currently only "magpie"

- target:

  name of a target dataset, see [`calcLandTarget`](calcLandTarget.md)
  for available target datasets

- harmonizationPeriod:

  Two integer values, before the first given year the target dataset is
  used, after the second given year the input dataset is used, in
  between harmonize between the two datasets

- harmonization:

  harmonization method, see [`toolGetHarmonizer`](toolGetHarmonizer.md)
  for available methods

## Value

harmonized wood harvest area data

## Author

Pascal Sauer

## Examples

``` r
if (FALSE) { # \dontrun{
  calcOutput("WoodHarvestAreaHarmonized", input = "magpie",
             target = "luh2mod", harmonizationPeriod = c(2015, 2050),
             harmonization = "fade")
} # }
```
