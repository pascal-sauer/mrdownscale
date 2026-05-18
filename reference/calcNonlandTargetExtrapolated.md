# calcNonlandTargetExtrapolated

Aggregated low resolution target data is extrapolated to the given years
using toolExtrapolate. To extrapolate wood harvest weight (bioh)
multiply wood harvest area already extrapolated by
calcLandTargetExtrapolated with the historical wood harvest rate in kg C
per Mha. Fertilizer is extrapolated and returned in kg ha-1 yr-1.

## Usage

``` r
calcNonlandTargetExtrapolated(input, target, harmonizationPeriod)
```

## Arguments

- input:

  name of an input dataset, currently only "magpie"

- target:

  name of a target dataset, see [`calcLandTarget`](calcLandTarget.md)
  for available target datasets

- harmonizationPeriod:

  Two integer values, will extrapolate to all years present in input
  data between harmonization start and end year

## Value

extrapolated nonland target data

## Author

Pascal Sauer

## Examples

``` r
if (FALSE) { # \dontrun{
  calcOutput("NonlandTargetExtrapolated", input = "magpie",
             target = "luh3", harmonizationPeriod = c(2025, 2050))
} # }
```
