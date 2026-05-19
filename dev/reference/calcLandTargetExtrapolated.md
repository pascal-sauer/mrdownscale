# calcLandTargetExtrapolated

Aggregated low resolution target data is extrapolated to the given years
using toolExtrapolate and normalized afterwards, so that the total sum
over all land types is unchanged. To account for the relationship
between wood harvest area and primary land (which is no longer primary
once it has been harvested) wood harvest area is calculated here even
though it is a nonland variable. The share of woody land that was
harvested in the historical period is calculated and then multiplied by
land (already extrapolated). Primary land is then converted to secondary
land, so that total reduction equals harvested area.

## Usage

``` r
calcLandTargetExtrapolated(input, target, harmonizationPeriod)
```

## Arguments

- input:

  name of an input dataset, options: "magpie", "witch"

- target:

  character, name of the target data set

- harmonizationPeriod:

  Two integer values, will extrapolate to all years present in input
  data between harmonization start and end year

## Value

extrapolated land target data, if calcOutput is called with
supplementary = TRUE and target is luh2mod wood harvest area is also
returned

## Author

Pascal Sauer
