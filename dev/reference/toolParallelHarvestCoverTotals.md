# toolParallelHarvestCoverTotals

For the absoluteChanges harmonizer, keep the global wood harvest area
totals of the two covers (forest: primf/secdf/pltns, other: primn/secdn)
parallel to the input (constant offset) after the harmonization year.

## Usage

``` r
toolParallelHarvestCoverTotals(
  xInput,
  xTarget,
  primHarv,
  landHarmonized,
  harmonizationYear
)
```

## Arguments

- xInput:

  recategorized wood harvest area input data, containing the
  harmonization year and the years after it

- xTarget:

  extrapolated wood harvest area target data, providing the anchor in
  the harmonization year

- primHarv:

  land-derived primary harvest area for the years after the
  harmonization year, with the items "primf" and "primn"

- landHarmonized:

  harmonized land data, providing the harvestable secondary land

- harmonizationYear:

  integer year the input changes are measured from

## Value

list with the secondary and plantation harvest area of the years after
the harmonization year (`secdHarv`) and the excess harvest area that
does not fit even in the other cover (`excessSecdHarvest`); the tool
checks its own allocation for regional and global parallelism of every
cover to the input

## Details

Processing steps, per year after the harmonization year:

- Regional cover target: anchor from the target data in the
  harmonization year plus the input change relative to that year.

- Set negative regional targets to zero. This raises their sum. To
  compensate, scale the regional targets \> 0 down, so the global cover
  total equals the raw global target. Regions set to zero get nothing.

- Set a negative global cover target to zero. The clipped loss is
  compensated for by the other cover: its regional targets are scaled by
  one common factor, so the global total over both covers stays input
  parallel. The scaled targets are the balanced targets, the resulting
  global total the balanced global target.

- Keep the primary harvest land-derived and distribute the remaining
  target over the clusters (weight: their harvestable secondary land),
  so the physical yearly cap for primary harvest holds by construction.

- Move harvest above the own cover's cap into the headroom of the other
  cover, region first and then globally. Only what fits nowhere becomes
  excess.

Output properties: In the ideal case, the harvest area of each cover
stays parallel to the input after the harmonization year, regionally and
globally, and so does the sum over both covers. When that is impossible,
parallelism is broken as follows:

1.  A negative regional target is set to zero and the targets \> 0 of
    the other regions are scaled down. The global cover total stays
    parallel, the parallelism of that region is lost.

2.  A cover whose global target is negative is set to zero, and the
    clipped loss is compensated for by the other cover. The global total
    of the cover set to zero is no longer parallel, the sum over both
    covers is.

3.  Where harvestable land is too little, harvest moves into the other
    cover's headroom, region first and then globally: first regional per
    cover parallelism is lost, then global per cover parallelism, the
    sum over both covers holds.

4.  The sum over both covers breaks where land-derived primary harvest
    alone outweighs a cover's balanced target, harvest exceeds anchor
    plus delta.

5.  Stranded cover targets of regions without any harvestable land for
    that cover fill the remaining room of all clusters
    (toolAllocateRegionTotals); what does not fit there is reported as
    excess directly, not shiftable.

Always true: the allocated secondary and plantation harvest stays within
the harvestable land of every cluster.

## Author

Pascal Sauer
