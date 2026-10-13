# plan: report deviations of the clamped (constantTotal = FALSE) absoluteChanges path

## Problem context

Goal: when plotting **global fertilizer totals** (`Tg yr-1`) of the `absoluteChanges`
harmonization against the unharmonized input
(`inst/extdata/plotting/plot-harmonization-nonland.R`, panel `fertilizer`), the harmonized
curve should follow the *raw* absolute-changes curve (= target of the harmonization year plus
the input's absolute changes, negatives untouched) as closely as possible while containing no
negative values.

Why the two curves differ (verified by reading the code, not yet quantified):

- Fertilizer is harmonized as **absolute Tg yr-1 per coarse crop category per cluster**
  (`R/calcNonlandHarmonized.R:55-60` builds `fertilizerTgTarget` via `toolFertilizerTg` and
  calls the harmonizer with `constantTotal = FALSE`). The negative handling for
  `constantTotal = FALSE` is a plain clamp at zero: `R/toolHarmonizeAbsoluteChanges.R:111-116`.
- Clamping **creates mass**: `sum(max(0, raw)) >= sum(raw)`, one-sided, growing with the
  distance to the harmonization year. So the harmonized global curve declines more slowly
  than raw (flattening, concave) — it sits *above* raw.
- Land does **not** have this problem, because for `constantTotal = TRUE` the negative
  handling is mass conserving (`toolHandleNegatives` group compensation,
  `R/toolHarmonizeAbsoluteChanges.R:66-110`) and it is already *reported* via
  `toolReportAreaDeviation(raw, out, groups)` at `R/toolHarmonizeAbsoluteChanges.R:108`.
  The clamped path has neither conservation nor reporting.
- Harmonization (land and nonland) runs at the **input/cluster resolution**, not at the high
  resolution grid (`calcLandHarmonized` uses `LandInputRecategorized`/`LandTargetLowRes`,
  `calcNonlandInputRecategorized` only attaches geometry). So "cells" here are clusters.
- Root cause of the negatives themselves: input fertilizer is regional MAgPIE Tg
  disaggregated by cropland (`R/calcNonlandInput.R:50-58`) while the target is LUH aggregated
  to clusters — cluster-level deltas mix two independent area trends, so clusters with
  `input_hy >> target_hy` clamp almost immediately.
- The plot's roundtrip is loss-free in one place and lossy in another:
  `toolFertilizerTg(rate, LandHarmonized)` inverts
  `toolFertilizerKgPerHa(out[, , "fertilizer"], landHarmonizedMha)`
  (`R/calcNonlandHarmonized.R:65`) with the same land object, so mass changes inside the
  harmonizer show up 1:1 in the plot. But the rate conversion separately loses absolute
  fertilizer where a category's harmonized cropland is 0 while its Tg is > 0 (`Inf` rates set
  to 0, `R/calcNonlandHarmonized.R:66-73`, documented by
  `tests/testthat/test-calcNonlandHarmonized.R:81`). That loss happens after the harmonizer and
  is **not** visible to the report being added here — keep it in mind when reading the numbers.

Missing measurement: unlike land, the clamped path reports nothing, so it is unknown how large
the clamp excess is and whether fertilizer dominates over `bioh` / `wood_harvest_area` /
`harvest_weight_type`. That measurement is this task's deliverable, and the precondition for
deciding how to close the gap in the plot.

## Task

Generalize and rename `toolReportAreaDeviation` into a unit-neutral report, and call it from
the clamped branch too, so the clamped nonland path reports how much value it added relative
to raw.

### 1. rename + generalize

`git mv R/toolReportAreaDeviation.R R/toolReportDeviation.R`, function becomes
`toolReportDeviation(xRaw, xOut, groups, unit = NULL)`:

- `unit = "Mha"` keeps the land column headers `=Mha/+Mha/-Mha`; `unit = NULL` means "values in
  the unit of each category", headers become `= / + / -` and the message wording must read
  accordingly (column suffix via `paste0(c("=", "+", "-"), unit)`).
- drop area-specific wording: "as % of the output area" → "as % of the output"; update the
  roxygen block, add `@param unit`. Keep everything else (sorting by gross deviation,
  decimal alignment, `%` of output accumulated over the post-harmonization years,
  `stopifnot(identical(getItems(xRaw, 3), getItems(xOut, 3)))`).
- new behaviour: skip a group whose positive, negative **and** net deviations are all below the
  existing `10^-5` tolerance, in addition to the existing empty-group `next`. This keeps the
  nonland run silent when nothing was clamped. It is an intentional behaviour change for land
  too (an untouched land group no longer prints an all-zero table).
- the tool stays internal (no `@export`) — the current one is not in `NAMESPACE` either, and
  tests reach internal functions through the package namespace. Do not edit `NAMESPACE` or
  `R/imports.R` (generated).

### 2. call it from the clamped branch

In `R/toolHarmonizeAbsoluteChanges.R`: rename the existing call at :108, and in the
`constantTotal = FALSE` branch (after `out[out < 0] <- 0`, :111-116) report with groups derived
automatically from the category prefix, because the nonland item set mixes units:

```r
toolReportDeviation(raw, out, groups = split(getItems(raw, 3), sub("\\..*$", "", getItems(raw, 3))))
```

Resulting groups: `wood_harvest_area`, `kgC_per_Mha`, `harvest_weight_type`, `fertilizer`.
`kgC_per_Mha` and `fertilizer` (Tg) are correct here — `calcNonlandHarmonized` hands the
harmonizer the substituted quantities, not the plotted ones (`R/calcNonlandHarmonized.R:27-60`).
Contract: `xRaw` holds the years after the harmonization year, `xOut` all years (true in both
branches); the tool slices `xOut[, yearsAfter, ]` internally.

### 3. tests (TDD: red first, then implement, then green)

- `git mv tests/testthat/test-toolReportAreaDeviation.R tests/testthat/test-toolReportDeviation.R`
  and rename the calls. The 4 existing tests assert only `group "cropland"` and numeric
  columns, so the wording change is inert for them.
- Move the `captureConditions` helper (`tests/testthat/test-toolReportAreaDeviation.R:2-16`) to
  `tests/testthat/helper-conditions.R` (testthat auto-loads `helper*.R`), since two files need
  it.
- Write the new report tests **before** the rename/implementation, referencing
  `toolReportDeviation` so they start red with "object not found":
  - `unit = NULL` produces `=`/`+`/`-` headers and the unit-neutral wording.
  - a group with only sub-tolerance deviations produces no message.
- `tests/testthat/test-toolHarmonizeAbsoluteChanges.R:296`
  (`"toolHarmonizeAbsoluteChanges works with nonland data"`) already exercises the clamp
  (`bioh.animal` raw `-1` in 2025 and `-2` in 2030 -> `0`):
  - assert the note now appears for `group "bioh"`, mean `+1.5` per timestep
    (`+3` accumulated over 2 timesteps) and group total `+8.57 %` (accumulated output
    `35` = `15 + 20` of `bioh.indu`, `0` of `bioh.animal`). Verify the exact formatted strings
    against real output when running instead of trusting these hand-computed digits; note
    `bioh.animal`'s own `%` is `Inf` (zero output).
  - add a second fixture that clamps nothing and assert **no** `bioh` note is emitted.
- Check the extra note does not break `tests/testthat/test-calcNonlandHarmonized.R`
  (`expect_no_warning` at :62, `expect_message(..., "Fertilizer application")` at :135).

### 4. verification

- single file: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-toolReportDeviation.R")'`
- `make format` then `make test` (all files must pass).
- real-data check (needs maintainer cache/data, not available in the sandbox): run
  `inst/extdata/plotting/plot-harmonization-nonland.R` and read the new `note` about deviations
  per nonland group.
- madrat cache: the report is a `tool*` function called inside the harmonizer body, and
  `toolGetHarmonizer.R:16-19` deliberately wraps harmonizers as `function(...) toolX(...)` so
  madrat records the dependency — changing these sources changes the cache keys of
  `LandHarmonized` / `NonlandHarmonized`, so cached outputs recompute automatically.

## Constraints

- never edit anything under `./man/` (generated). `man/toolReportAreaDeviation.Rd` becomes
  orphaned by the rename; leave it to the maintainer (`make docs`).
- do not run `make install`, `make docs`, `make build`.
- style: `if`/`for`/`while` need `{` and a newline (single-line ternary style `if (cond) 1 else 2`
  is the only exception); one function per file, filename == function name; no comments unless
  they carry non-obvious information.
