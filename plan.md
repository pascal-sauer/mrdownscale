# Plan: implement `harmonization = "absoluteChanges"` in calcLandHarmonized

## Context

`calcLandHarmonized` currently supports three harmonization methods (`fade`,
`fadeForest`, `offset`), selected by name via `toolGetHarmonizer()` and applied as
`harmonizer(xInput, xTarget, harmonizationPeriod)`. We need a fourth method,
`absoluteChanges`, that applies the *absolute changes* of the input data (relative
to a single harmonization year) on top of the historical target value at that year:

```
out[t] = target[harmonizationYear] + (input[t] - input[harmonizationYear])   for t >= harmonizationYear
out[t] = target[t]                                                           for t <  harmonizationYear
```

Example: harmonizationYear = 2020, target secdf(2020) = 10 Mha, input secdf grows
by 2 Mha from 2020→2025 ⇒ out secdf(2025) = 12 Mha.

This is mathematically the aneris/`mip` "offset" method **without convergence** (a
permanent offset). Because input total area is constant over time, the total area is
automatically conserved; only *individual* categories can go negative.

Key difference from existing methods: `harmonizationPeriod` is a **single year**
(here called `harmonizationYear`) that must exist in the target data, not a
`c(start, end)` pair. This means the existing pipeline plumbing that relies on
`hp[2]` (e.g. `LandTargetExtrapolated`, and the post-harmonization "returns input
after hp[2]" checks) does not apply and must be branched around.

### Decisions (confirmed with user)
- **Negative categories**: clamp to 0, then renormalize each region/timestep so the
  total area is preserved; additionally emit a **warning** when any clamping happens.
  Also run `toolReplaceExpansion` for `primf`/`primn` as the other harmonizers do.
- **Scope**: make it work end-to-end for `calcLandHarmonized` only. Do **not** touch
  `calcNonlandHarmonized` / `calcWoodHarvestAreaHarmonized` (they use `hp[2]`).

## Reference material (existing patterns to reuse)

- `R/toolGetHarmonizer.R` — registry `list(offset=..., fade=..., fadeForest=...)`;
  wrap the new function as `absoluteChanges = function(...) toolHarmonizeAbsoluteChanges(...)`
  so madrat still detects it as a dependency (comment in that file explains why).
- `R/toolHarmonizeFade.R` — canonical harmonizer structure: `stopifnot` validation of
  years/items, `xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]`
  alignment, `mbind` of time slices, and `toolReplaceExpansion(out, "primf", "secdf", ...)`
  / `(out, "primn", "secdn", ...)`.
- `R/calcLandTargetExtrapolatedCore.R` normalization idiom for conserving total area:
  `x <- x * targetArea / dimSums(x, dim = 3)` (used here after clamping).
- `R/calcLandTargetLowRes.R` — `calcOutput("LandTargetLowRes", input, target,
  endOfHistory = <year>)` returns low-res target land through a given year; use this
  to source historical target for the single-year path (avoids `LandTargetExtrapolated`'s
  `hp[2]` dependency).
- `toolStatusMessage("warn", ...)` / `toolExpectTrue(..., falseStatus = "warn")` — the
  codebase's warning idioms (see `toolEqualizeArea.R`, `toolPrimExpansionCheck.R`).
- `CONTRIBUTING.md`: define aux functions inside the `calc`/`tool` body or as a `tool*`
  function so they enter madrat's cache key. The new harmonizer is a `tool*` function,
  so this is satisfied.

## Implementation steps (strict TDD — for each code step: write test, run it, confirm
it FAILS, implement, run again, confirm it PASSES; do not advance until green)

### Step 1 — New harmonizer `toolHarmonizeAbsoluteChanges` (unit-tested in isolation)

New file `R/toolHarmonizeAbsoluteChanges.R`. Signature mirrors the others:
`toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod, level = 3)`.

Behavior:
1. `y <- harmonizationPeriod`; `stopifnot(length(y) == 1, round(y) == y)`.
2. Validate: `y %in% getYears(xTarget, TRUE)`, `y %in% getYears(xInput, TRUE)`,
   region/item set-equality between input and target (same `stopifnot` style as fade).
3. Align: `xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]`.
4. `futureYears <- inputYears[inputYears > y]`.
5. `delta <- xInput[, futureYears, ] - setYears(xInput[, y, ], NULL)` (broadcast the
   single harmonization-year slice across future years).
6. `changed <- setYears(xTarget[, y, ], NULL) + delta` → set years to `futureYears`.
7. `out <- mbind(xTarget[, targetYears <= y, ], changed)`.
8. `out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100, level = level)`
   and same for `primn`/`secdn` (guard with the `%in% getItems` checks as fadeForest does).
9. Negative handling: if `any(out < 0)`, emit `toolStatusMessage("warn", ...)` (or
   `toolExpectTrue(min(out) >= 0, ..., falseStatus = "warn")`), clamp `out[out < 0] <- 0`,
   then renormalize per region/timestep to restore the total:
   `out <- out * setYears(dimSums(xTarget[, y, ], 3), NULL) / dimSums(out, 3)` guarding
   against divide-by-zero (`out[is.na(out)] <- 0`), matching the Core normalization idiom.
10. Return `out` (magpie), sets/order consistent with inputs.

Test `tests/testthat/test-toolHarmonizeAbsoluteChanges.R` (mirror the fixture style of
`test-toolReplaceExpansion.R`, build objects with `new.magpie`/`getItems`/`getYears`):
- **core formula**: target years {2000,2010,2020}, input years {2020,2025,2030},
  `harmonizationPeriod = 2020`; assert output years = {2000,2010,2020,2025,2030},
  `out[,2020,] == target[,2020,]`, `out[,2025,] == target[,2020,] + (input[,2025,] - input[,2020,])`.
- **total conserved**: `dimSums(out,3)` equals target total at 2020 for every year.
- **pre-year passthrough**: `out[, year < 2020, ] == target[, year < 2020, ]`.
- **negative clamp + warning**: craft an input where one category's decrease exceeds
  the target baseline; assert `expect_warning(...)`, `all(out >= 0)`, and total still
  conserved.

### Step 2 — Register the method in `toolGetHarmonizer`

Edit `R/toolGetHarmonizer.R`: add
`absoluteChanges = function(...) toolHarmonizeAbsoluteChanges(...)` to the list; update
the roxygen `@param`/`@seealso` to mention it.

Test: extend/add `tests/testthat/test-toolGetHarmonizer.R` asserting
`is.function(toolGetHarmonizer("absoluteChanges"))` and that an unknown name still errors.

### Step 3 — Branch `calcLandHarmonized` for the single-year method

Edit `R/calcLandHarmonized.R`:
- Introduce `hp1 <- harmonizationPeriod[1]` (works for length-1 and length-2).
- Source the target conditionally:
  - `absoluteChanges`: `stopifnot(length(harmonizationPeriod) == 1)`; get
    `xTarget <- calcOutput("LandTargetLowRes", input = input, target = target,
    endOfHistory = harmonizationPeriod, aggregate = FALSE)` and assert
    `harmonizationPeriod %in% getYears(xTarget, TRUE)`.
  - otherwise: keep the existing `LandTargetExtrapolated` call.
- The area checks (lines 24–31) and `toolEqualizeArea(xInput, xTarget[, hp1, ])` stay as
  is (use `hp1`).
- The "Returning reference data before harmonization period" check (out[<=hp1] ==
  target[<=hp1]) stays — it holds for absoluteChanges.
- **Guard the post-harmonization input-equality block (current lines 55–68) with
  `if (harmonization != "absoluteChanges")`** — those checks assume `out == input`
  after `hp[2]`, which is false (and `hp[2]` does not exist) for absoluteChanges. The
  method's correctness is instead covered by the harmonizer's own unit tests plus the
  retained common checks (`all(out >= 0)`, totals constant, `toolPrimExpansionCheck`).
- Update the roxygen `@param harmonizationPeriod` note to document that for
  `absoluteChanges` it is a single year existing in the target.

Because a full `calcLandHarmonized` unit test needs madrat data fixtures (the two
existing tests are tool-level only), this step is validated by the Step 4 end-to-end run
and its internal `toolExpect*` checks rather than a new isolated test.

### Step 4 — End-to-end verification (definition of done)

In a **fresh** R session:

```r
pkgload::load_all()
a <- calcOutput("LandHarmonized", input = "magpie", target = "luh3",
                harmonizationPeriod = 2020, harmonization = "absoluteChanges",
                aggregate = FALSE)
```

Must complete without error. If the harmonized object is cached, force re-execution with
`madrat::setConfig(ignorecache = c("calcLandHarmonized", "calcLandTargetLowRes"))`.
Spot-check: `a[,2020,]` matches the target at 2020, and post-2020 changes equal the
input's changes (up to any clamping), and `dimSums(a,3)` is constant over time.

## Files

- New: `R/toolHarmonizeAbsoluteChanges.R`
- New: `tests/testthat/test-toolHarmonizeAbsoluteChanges.R`
- Edit: `R/toolGetHarmonizer.R` (registry + roxygen)
- New/Edit: `tests/testthat/test-toolGetHarmonizer.R`
- Edit: `R/calcLandHarmonized.R` (branch target sourcing + guard post-harmonization checks + roxygen)

## Verification summary

1. `Rscript -e 'devtools::test_file("tests/testthat/test-toolHarmonizeAbsoluteChanges.R")'` — green.
2. `Rscript -e 'devtools::test_file("tests/testthat/test-toolGetHarmonizer.R")'` — green.
3. `Rscript -e 'devtools::test()'` — no regressions.
4. Fresh-session `calcOutput("LandHarmonized", ..., harmonization = "absoluteChanges")` — succeeds.

Do not run `lucode2::buildLibrary()` / bump versions unless the user asks.

