# toolParallelCoverTargets

The cover target is the absolute cover per region and year: the value
from the target data in the harmonization year plus the change in the
input since that year.

## Usage

``` r
toolParallelCoverTargets(inputCover, targetCover, harmonizationYear)
```

## Arguments

- inputCover:

  named list of region-level per-year cover totals, the input including
  the harmonization year and the years after it

- targetCover:

  named list of region-level per-year cover totals, providing the anchor
  in the harmonization year

- harmonizationYear:

  integer year the input changes are measured from

## Value

list with the balanced regional cover targets per cover (`coverTargets`)
and the balanced global cover total (`coverGlobalBalanced`)

## Details

Per cover, the tool sets negative regional targets to zero. This raises
their sum. To compensate, it scales the targets \> 0 of the same cover
in the other regions down, so the global total equals the raw
(unmodified) global target. Covers do not interact here.

The tool also sets a negative global cover target to zero. It
compensates for the clipped loss by scaling the regional targets of all
covers with one common factor, so the global total over all covers is
met (zero if it is negative). Returns the balanced regional cover
targets and the balanced global cover total.

## Author

Pascal Sauer
