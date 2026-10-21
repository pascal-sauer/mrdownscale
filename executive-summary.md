# Executive summary

## Regional compensation of lost fertilizer tonnage

The `absoluteChanges` compensation for fertilizer tonnage lost to zero clamping
and vanished cropland used one global factor per crop group and year, so losses
in one region inflated rates everywhere else. It now rescales per region, with
the global factor only as fallback in crop group x year tables where a region
lost all its cropland and can no longer hold its expected tonnes (reported via
a warn check).

Production verify (magpie/luh3, harmonization 2020, measures in
`/p/tmp/R_tmp/opencode/fert-regional/`): the warn fires once for `NEU
fertilizer.c3nfx 2040` (0.0033 Tg unholdable), the global parallelism check and
`toolCheckFertilizer` stay green, and 753/900 region totals shift back toward
their own region (total 135 Tg, max 2.14 Tg) while the plausibility cap no
longer needs to move tonnes across regions (57.5 Tg -> ~0).
