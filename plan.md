# Plan: handle stranded bioh tonnage in cells with zero harmonized harvest area

Branch: `fix-bioh-normalization` (off `fix-wood-harvest` @ ad913bb). Do not implement before reading this fully.

## Problem

6 ScenarioMIP downscaling runs (`/p/tmp/pascalfu/scenariomip-downscaling`, `Rscript downscale.R <scenario>`, see `slurm.sh`). 4 of 6 halt:

```
Error in mrdownscale:::calcNonlandHarmonized(...) :
  biohNormalization < Inf are not all TRUE
Calls: retrieveData ... -> stopifnot
```

Failing: SSP3-H, SSP2-M, SSP2-ML, SSP2-LN (logs `slurm-2520139/40/41/43.log`, same in `old-slurm-logs/`). Passing: SSP2-L, SSP1-VL2. Failing runs produce no `multiple-management`/`multiple-transitions` output.

## Root cause

`R/calcNonlandHarmonized.R:107-129` adapts harmonized bioh tonnage to the harmonized wood harvest area:

- `kgCPerMhaHarmonized` (harmonized kg C / Mha rates) is floored strictly > 0 (`:109-114`).
- `biohCalculated <- kgCPerMhaHarmonized * collapseDim(harvestArea)` (`:115`) is therefore 0 iff the cell-year's **total** harmonized harvest area is 0.
- `biohNormalization <- dimSums(out[, , "bioh"], 3) / dimSums(biohCalculated, 3)` (`:119`) rescales per cell-year so the harmonized bioh total is kept (`stopifnot` `:121`, soft check `:125-127`). 0/0 -> 0 handled (`:120`), but **positive/0 -> Inf** hits `stopifnot(biohNormalization < Inf)` (`:121`).

`calcWoodHarvestAreaHarmonized` drives harvest area to 0 per cell via physical limits: `primHarv = pmax(prim land reduction / timestep, 0)` (`R/calcWoodHarvestAreaHarmonized.R:40-43`, 0 in afforesting cells) and secondary harvest capped at `toolMaxHarvestPerYear` = previous-timestep secondary/plantation area (`R/calcWoodHarvestAreaHarmonized.R:76-89`, `R/toolMaxHarvestPerYear.R:27`, 0 in cells without secondary land). Harmonized bioh tonnage stays > 0 because `toolHarmonizeAbsoluteChanges` anchors it to the LUH3 value of 2020 and only clamps to 0 (`R/toolHarmonizeAbsoluteChanges.R:115-116`). In steep-decline scenarios MAgPIE collapses a cell's harvest *area* faster than its *tonnage* (intensification), so both occur in the same cell (log evidence: `slurm-2520143.log:96-126`, deviation tables show clamp deviations only positive, e.g. `bioh.primn +197%`). Such cells cannot satisfy both invariants (tonnage kept AND tonnage = rate x area) - the tonnage is stranded and must move.

## Fix (agreed: option 1 - redistribute stranded tonnage, no plausibility cap)

Analogous to fertilizer handling (`R/calcNonlandHarmonized.R:75-102`, `R/toolCapFertilizerKgPerHa.R`): rates x land cannot host tonnes where land vanished -> tonnes move to cells that still have land, region first, global pool fallback. Here: rates have no plausibility cap (harvest yield varies too much to pick one - deliberate), so tonnage simply moves to cells with harvest area > 0.

Design:

1. New internal tool `R/toolRedistributeBioh.R` (not exported in `NAMESPACE`, roxygen docs per `CODE_STYLE.md`, modeled on `toolCapFertilizerKgPerHa.R`):
   - `toolRedistributeBioh(biohTonnage, biohCalculated)` where both share full item sets `bioh.*` (6 subcategories), cell x year.
   - receivers = cells with `collapseDim(biohCalculated, 3) > 0`, weight = `biohCalculated` per subcategory per cell-year.
   - per subcategory and year: stranded tonnage from cells with `collapseDim(biohCalculated) == 0` is added to receivers of the **same region** proportional to their weight; if a region has no receivers, use the **global pool**; if nowhere has receivers, keep 0 and report via `toolExpectTrue` (soft, like the infeasible fertilizer fallback `R/toolCapFertilizerKgPerHa.R:33-35`).
   - single pass suffices (receivers only receive, no overshoot); tonnage totals exactly preserved except in the reported no-receiver-anywhere case.
2. In `R/calcNonlandHarmonized.R`, apply it to the bioh items of `out` **before** the normalization block (before `:119`):
   - stranded mask: `collapseDim(harvestArea) == 0` (after `:105`) intersected with positive `out[, , "bioh"]` totals.
   - `out[, , "bioh"] <- toolRedistributeBioh(out[, , "bioh"], biohCalculated)`; stranded cells become exactly 0 there, so at `:119` numerator and denominator are 0 together and `:120`'s `is.nan -> 0` keeps working. Normalization `:119-127` otherwise unchanged - per-cell total-keeping check `:125` still passes by construction.
   - region totals of bioh are preserved, so the global bioh trajectory stays parallel to the MAgPIE input.
3. Docstring of `calcNonlandHarmonized` (lines 6-9): mention that tonnage of cells without harvest area is redistributed within region/global.

## Reproduction (write failing test first, then fix, then green)

Add to `tests/testthat/test-calcNonlandHarmonized.R`, reusing the mock pattern of the first test (`local_mocked_bindings(calcOutput = ...)` with `xInput`/`xTarget`/`landTarget`/`landHarmonized`/`harvestArea`):

- fixture identical to "calcNonlandHarmonized works with harmonization = absoluteChanges", except `harvestArea["reg.one", c(2025, 2030), ] <- 0` (cell keeps positive bioh tonnage after 2020, zero harvest area; `reg.two` keeps area 1 everywhere).
- pre-fix it reproduces the crash: `expect_error(calcNonlandHarmonized(...), "biohNormalization < Inf")`.
- post-fix assertions:
  - no error, no "Check failed" condition (use `captureConditions` helper from `tests/testthat/helper-conditions.R`).
  - `x["reg.one", c(2025, 2030), paste0("bioh.", c("primf", "secmf"))] == 0`.
  - reg.two's bioh gained exactly reg.one's stranded tonnage per subcategory (region totals equal the pre-fix harmonizer totals; hand-compute expected values from the fixture like the existing tests do).
  - `all(is.finite(x))`, `all(x >= 0)`, fertilizer/other items unchanged vs the first test.
- keep a second case? Only needed if trivial: region without any receiver falls back to global / reports (a tool-level test of `toolRedistributeBioh` mirroring the two `toolCapFertilizerKgPerHa` tests at the end of the file is enough).

Real-data repro for cluster validation (write as `/p/tmp/pascalfu/scenariomip-downscaling/repro-bioh.R`, NOT inside the package repo; needs the cluster env - `/p/projects/scenarioMIP` is not listable from the login node, run via the same hook as `slurm.sh`):

```r
library(mrdownscale)
setConfig(outputfolder = "./madrat-outputfolder", maxLengthLogMessage = 500)
redirectSource("MagpieFulldataGdx",
               target = file.path("/p/projects/scenarioMIP/main_v8p5/remind/magpie/output",
                                  "C_SMIPv08-VLHO-SSP2-EcPrice1300-def-mag-5",
                                  c("fulldata.gdx", "clustermap_rev4.119_c200_67420_h12.rds")),
               linkOthers = FALSE)
x <- calcOutput("NonlandHarmonized", input = "magpie", target = "luh3",
                harmonizationPeriod = c(2020, 2020), harmonization = "absoluteChanges",
                aggregate = FALSE)
stopifnot(all(is.finite(x$x)))
message("max: ", max(x$x))
```

With existing madrat caches it errors within seconds (log `slurm-2520143.log:86-136`); a cache miss re-parses the ~600 MB gdx (minutes). Diagnostic option: count/sum stranded cells via `collapseDim(calcOutput("WoodHarvestAreaHarmonized", input = "magpie", target = "luh3", harmonizationPeriod = c(2020, 2020), harmonization = "absoluteChanges", aggregate = FALSE)) == 0`.

## Verification

1. `make test` in the package dir (devtools::test). Note: local renv is out-of-sync and its madrat 3.44.0 lacks `toolExpectTrue`/`toolExpectLessDiff` - use the same R environment the slurm hook (`RSCRIPT_SLURM_HOOK`) uses if the local one fails to load the package.
2. Re-run all 6 scenarios via `slurm.sh` on the cluster: the 4 failing ones must now finish; the 2 previously passing (SSP2-L, SSP1-VL2) must produce unchanged output (no stranded cells there -> redistribution is a no-op; spot-check via the consistencyCheck.log and file checksums).
3. Watch the `bioh` deviation report (`toolReportDeviation`) - magnitudes must not change (redistribution happens after the harmonizer).
4. Final step: write a short, concise `executive-summary.md` in the worktree folder:
   - what was done, why, how (few sentences each, no filler);
   - how it improved the situation, as a before vs after comparison: before, 4 of 6
     ScenarioMIP runs (SSP3-H, SSP2-M, SSP2-ML, SSP2-LN) crashed with
     `biohNormalization < Inf are not all TRUE` and produced no
     `multiple-management`/`multiple-transitions` output; after, all 6 runs finish,
     the 2 previously passing ones unchanged, bioh region/global totals preserved.

## Rejected alternatives (do not re-litigate)

- `biohNormalization[!is.finite(.)] <- 0` + relaxed check: silently drops harvest carbon vs MAgPIE input in 4/6 scenarios; downstream carbon accounting reconciles bioh.
- Derive bioh purely from kgC/Mha rates x area (drop renormalization): changes global bioh for all scenarios.
- Plausibility cap on kgC/Mha for receivers: no defensible threshold, deferred as follow-up.
- Flooring harvest area upstream: zeros come from physical land limits, would break `toolCheckWoodHarvestArea`/excess-harvest consistency.
