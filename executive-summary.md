# Executive summary: harmonized fertilizer tonnages projected onto the harmonized land

Spatial vocabulary: *cells* are the 0.25° target grid cells (e.g. `-179p75.-82p75`),
*clusters* are the MAgPIE spatial units and dim-1 items of the harmonized data
(e.g. `JPN.58`), *regions* are the MAgPIE region codes (`JPN`, `CHA`, ...).

## Problem

The production run (`input = "magpie"`, `target = "luh3"`, `harmonizationPeriod =
c(2020, 2020)`, `harmonization = "absoluteChanges"`) failed `toolCheckFertilizer`:
192 cluster x year x crop fertilizer values exceeded the plausibility threshold of
1200 kg ha-1 yr-1, with a global maximum of 1.3e10. Root cause: fertilizer was
harmonized in tonnage space per cluster (target anchor plus MAgPIE input delta), but
MAgPIE intensifies crops that the harmonized (LUH3-based) land says are vanishing —
tonnes grew on clusters whose harmonized crop-group area collapsed, and were then
divided by that collapsed area.

## What was done

In the `absoluteChanges` branch of `calcNonlandHarmonized` the harmonized land is
now treated as truth for fertilizer tonnages:

- `toolProjectFertilizerTg` (new): the input tonnage change of each region since
  the harmonization year is redistributed within the region onto its clusters
  proportional to the harmonized crop-group areas of each year, on top of the
  unchanged per-cluster target anchor. Region sums are preserved, and new tonnes
  only reach clusters where the harmonized land has that crop. The existing
  global parallelism rescale then runs unchanged.
- `toolCapFertilizerKgPerHa` (new): one non-iterative residual pass caps rates at
  1200 kg ha-1 yr-1 (the `toolCheckFertilizer` threshold), moving each cluster's
  excess into the remaining headroom of its region's clusters for that crop group
  and year (global pool as fallback, input returned unchanged if infeasible so
  totals stay intact and the check reports it). This subsumes the abandoned
  `fert-area-floor` area-floor experiment — residue clusters are clamped and their
  tonnes relocated instead of zeroed.

TDD: red tests first (collapse fixture, residue fixture, region sum preservation,
cap redistribution/fallback/infeasible cases), then implementation; `make format`
and `make test` green (212 tests).

## Results (production cache, baseline recomputed from b5775e1)

| Measure | before | after |
| --- | --- | --- |
| clusters > 1200 kg ha-1 yr-1 | 192 | 0 |
| global max rate (kg ha-1 yr-1) | 1.30e10 | 1200 |
| toolCheckFertilizer | fails (max: 1.3e+10) | passes (max: 1200) |
| global parallelism check | passes | passes (deviation 4.3e-14) |
| pre-2020 output vs target | identical | identical (max diff 0) |
| post-hp total tonnage (Tg) | 1269.8 | 1269.8 |

Plausibility: 262.7 Tg (0.5 x sum |delta|, 21% of the post-2020 total) relocated
between clusters relative to before — that is the projection moving intensification
onto existing cropland plus the cap moving anchor residues off collapsed areas.
Region x crop-group mean rates shift by a median of 4.9 kg ha-1 yr-1; the largest
shifts (e.g. JPN c4ann, 481 -> 110) are negative regional deltas now applied
regionally instead of being kept by zero-clamped clusters. Rescale factors span
0.50-1.50 (previously 0.39-1.22), same mechanism and magnitude. Measurement
artifacts: `/p/tmp/R_tmp/opencode/fert-compare/` (`fert_before.rds`,
`fert_proj.rds`, `measure4.R`/`measure4.log`).
