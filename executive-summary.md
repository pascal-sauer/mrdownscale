# Executive summary: stranded bioh in clusters with zero harmonized harvest area

## What

New internal tool `R/toolRedistributeBioh.R`, applied in `calcNonlandHarmonized()` to
the harmonized bioh before its normalization step. Covered by two tool-level tests
and one `calcNonlandHarmonized()` test reproducing the production failure.

## Why

4 of 6 ScenarioMIP downscaling runs (SSP3-H, SSP2-M, SSP2-ML, SSP2-LN) crashed with
`biohNormalization < Inf are not all TRUE` and produced no
`multiple-management`/`multiple-transitions` output. In steep-decline scenarios
MAgPIE collapses a cluster's harvest *area* (physical limits: no land reduction to
spread, no secondary land to harvest) faster than its *bioh*, which the
`absoluteChanges` harmonizer anchors to the LUH3 value of 2020 and only clamps at
zero. Such clusters have positive harmonized bioh but zero recalculated bioh (kg C per
Mha rates times harvest area), so the normalization became positive / 0 = Inf.
They cannot keep their bioh and satisfy bioh = rates x area at the same time.

## How

Before the normalization, bioh on clusters without harvest area is moved to the clusters
that still have some: per subcategory and year to the receivers of the same region
proportional to their recalculated bioh, falling back to the global pool when a
region has no receivers. Bioh that cannot be moved anywhere (no receivers at all
for a subcategory and year) is dropped and reported as a warning, so the tool's
invariant holds: returned bioh is zero wherever recalculated bioh is zero (the
converse does not hold, zero-clamped changes legitimately leave bioh zero on clusters
with harvest area). Every redistribution is logged: absolute and relative amount
and affected cluster years. Receivers only receive, so one pass suffices and no one
overshoots; region and global bioh totals are preserved unless bioh had to be
dropped. Stranded clusters then hold exactly zero, so normalization sees 0/0
-> 0 as before. Like the analogous fertilizer handling
(`toolCapFertilizerKgPerHa`), but without a plausibility cap on the rates, which
harvest yield variation does not justify.

## Before vs after

- Before: 4 of 6 runs halted inside `calcNonlandHarmonized()`, no downstream output.
- After: the crash is reproduced by the new test pre-fix and gone post-fix; all
  package tests pass (246 PASS, lint clean); stranded clusters end at zero, receiving
  clusters gain exactly the stranded bioh, global totals equal the harmonizer totals;
  runs without stranded clusters (SSP2-L, SSP1-VL2) are a no-op, as are the bioh
  deviation reports, since redistribution happens after the harmonizer.
- Pending: re-run of all 6 scenarios via `slurm.sh` (the 4 failing ones must
  finish, the 2 passing ones must produce checksum-identical output);
  `repro-bioh.R` reproduces the SSP2-LN crash on real data within seconds and
  doubles as the post-fix check (`stopifnot(all(is.finite(x$x)))`).
