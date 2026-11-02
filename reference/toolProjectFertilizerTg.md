# toolProjectFertilizerTg

Projects fertilizer tonnages onto the harmonized land: for every year
after the harmonization year, the change of each region's tonnage total
relative to the harmonization year is distributed onto the clusters of
the region proportional to their harmonized crop group area of that
year, on top of the per-cluster target tonnage of the harmonization
year. Region totals are preserved, and tonnes only grow where the
harmonized land has cropland of the respective crop group.

## Usage

``` r
toolProjectFertilizerTg(
  fertilizerTgTarget,
  fertilizerTgInput,
  landMha,
  harmonizationYear
)
```

## Arguments

- fertilizerTgTarget:

  target fertilizer in Tg yr-1, provides the per-cluster tonnage anchor
  for the harmonization year

- fertilizerTgInput:

  input fertilizer in Tg yr-1, containing the harmonization year and the
  years after it

- landMha:

  harmonized land in Mha

- harmonizationYear:

  integer year the input changes are measured from

## Value

projected fertilizer in Tg yr-1 for the years after the harmonization
year, clamped at zero, in the item order of fertilizerTgTarget

## Author

Pascal Sauer
