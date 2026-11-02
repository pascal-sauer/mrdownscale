# toolCapFertilizerKgPerHa

Enforce a plausibility threshold on fertilizer rates by redistributing
the tonnage of every cell above it onto the remaining headroom of the
cells of the same region, crop group and year, or of the global pool for
those crop group and year tables in which a region cannot absorb its own
excess. Receivers only receive a share proportional to their own
headroom, so a single pass suffices and no one overshoots. Tables in
which not even the global pool has enough headroom are left unchanged,
so that the tonnage totals stay intact and toolCheckFertilizer reports
the exceedance.

## Usage

``` r
toolCapFertilizerKgPerHa(fertilizerKgPerHa, landMha, threshold = 1200)
```

## Arguments

- fertilizerKgPerHa:

  fertilizer rates in kg ha-1 yr-1

- landMha:

  land in Mha, used to convert to and from tonnages

- threshold:

  maximum plausible rate in kg ha-1 yr-1

## Value

fertilizer rates in kg ha-1 yr-1, capped at the threshold

## Author

Pascal Sauer
