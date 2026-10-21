# Plan: compensate lost fertilizer tonnage regionally instead of with one global scaling factor

## Context

This branch (`fert-land-projection`) already projects the regional input deltas
onto the harmonized land (`toolProjectFertilizerTg`) and caps anchor residues
(`toolCapFertilizerKgPerHa`) — see the superseded plan of that work, its
measurements in `/p/tmp/R_tmp/opencode/fert-compare/` and `executive-summary.md`.

Remaining gap: the compensation step in the `absoluteChanges` block of
`R/calcNonlandHarmonized.R:75-102`. Where the zero clamping of the harmonizer
and the `Inf -> 0` replacement (cells whose cropland vanished, i.e. the
"compensating negative fertilizer" step) drop tonnage, the lost tonnes are
restored by **one global scaling factor per crop group x year**:

```r
factor <- expectedGlobalTg / actualGlobalTg   # R/calcNonlandHarmonized.R:92
fertilizer[, yearsAfterHy, ] <- fertilizer[, yearsAfterHy, ] * factor
```

So a loss in region A inflates the fertilizer rates of regions B and C, even
though region A's own remaining cells could absorb it. This violates MAgPIE's
between-region structure (the very reason the projection in
`toolProjectFertilizerTg` is regional). The residual cap right below already
does it better: redistribute **within region**, global pool only as fallback
(`R/toolCapFertilizerKgPerHa.R:26-42`).

Goal: the rescale must follow the same pattern — per region x crop x year, so
lost tonnes are compensated by the other clusters of the SAME region; the
global parallelism check (`R/calcNonlandHarmonized.R:158-161`) stays exact as a
consequence.

## Design

Replace the global rescale in the second `absoluteChanges` block of
`calcNonlandHarmonized` (lines 75-94, everything from `inputGlobal` down to the
`fertilizer[, yearsAfterHy, ] <- ... * factor`, before the cap call) with:

```r
regionTgInput <- dimSums(fertilizerTgInput, 1.2)
expectedRegionTg <- setYears(dimSums(fertilizerTgTarget[, hp[1], ], 1.2), NULL) +
  (regionTgInput[, yearsAfterHy, ] - setYears(regionTgInput[, hp[1], ], NULL))
expectedRegionTg[expectedRegionTg < 0] <- 0
# the global check below stays exact: it compares against the regional
# expectations summed up
expectedGlobalTg <- dimSums(expectedRegionTg, 1)

# where the zero clamping of the harmonizer removed all cropland of a cell,
# the fertilizer rate is undefined, so set it to zero like toolFertilizerKgPerHa
# does for cells without fertilizer
fertilizer <- out[, , "fertilizer"]
fertilizer[fertilizer == Inf] <- 0

# rescale per region, crop group and year so that every region's total equals
# its target total plus its own input delta, which redistributes fertilizer
# lost in cells without cropland to the cells of the SAME region that still
# have some
actualRegionTg <- dimSums(toolFertilizerTg(fertilizer[, yearsAfterHy, ], landHarmonizedMha), 1.2)
if (!toolExpectTrue(!any(expectedRegionTg[] > 0 & actualRegionTg[] == 0),
                    "Regions with expected fertilizer tonnage have cropland left to hold it",
                    falseStatus = "warn")) {
  # at least one region lost all cropland of a crop group but still expects
  # tonnes, it can never hold them, so fall back to the global pool, whose
  # total is what the parallelism check constrains
  actualGlobalTg <- dimSums(actualRegionTg, 1)
  factor <- expectedGlobalTg / actualGlobalTg
  factor[actualGlobalTg == 0] <- 1
} else {
  factor <- expectedRegionTg / actualRegionTg
  factor[actualRegionTg == 0] <- 1
}
fertilizer[, yearsAfterHy, ] <- fertilizer[, yearsAfterHy, ] * factor
```

Notes:

- `fertilizer * factor` with a region-level factor broadcasts onto the cells
  via magclass region-set matching — established precedent:
  `R/toolCapFertilizerKgPerHa.R:37-43` multiplies a `dimSums(..., 1.2)` factor
  with cell-level tonnage.
- `toolExpectTrue` (mstools) invisibly returns the check result and reports a
  false check as `"Check failed: <description>"` via `toolStatusMessage` with
  the given `falseStatus` — house style, see e.g. `R/calcLandTarget.R:105-106`.
  Maintainer decision: use this, NOT a plain `warning()`.
- The per-region clamp (`< 0 -> 0`) replaces the old global clamp of
  `expectedGlobalTg`; `expectedGlobalTg <- dimSums(expectedRegionTg, 1)` keeps
  the existing global parallelism check (lines 158-161, unchanged code) exact
  in BOTH branches (fallback factor multiplies global actual to exactly that
  total).
- `actualRegionTg == 0 & expectedRegionTg == 0` (region never had any) must
  give factor 1, hence the `factor[...] <- 1` guards.

Strengthened check — inside the final `absoluteChanges` check block
(`R/calcNonlandHarmonized.R:158-161`), BEFORE the existing global one:

```r
toolExpectLessDiff(dimSums(toolFertilizerTg(out[, yearsAfterHy, "fertilizer"], landHarmonizedMha), 1.2),
                   expectedRegionTg, 10^-5,
                   "Regional fertilizer totals are parallel to input after the harmonization year")
```

In the fallback case this per-region check fails by design (that is the whole
point of the fallback) and reports a second `"Check failed"` note — accepted
(maintainer decision 2a: always report, no skipping).

Update the roxygen description of `calcNonlandHarmonized`
(`R/calcNonlandHarmonized.R:10-18`): the rates are rescaled **per region** to
keep each region's total parallel to the input, with a global fallback if a
region lost all cropland of a crop group. Do NOT touch `man/` (maintainer
regenerates).

## Verified fixtures

`/p/tmp/R_tmp/opencode/fert-regional/draft.R` contains both fixtures already
run against the current code (measured values below are real). All tonnages
use `Tg = kg ha-1 * Mha * 10^-3`.

**Fixture 1 — regional compensation** (cells `regone.1`, `regone.2`,
`regtwo.1`; harmonization year 2020; `c3ann` only):

- `landHarmonized`: 100 Mha everywhere, `regone.2` -> 0 in 2025.
- target rate 100 kg ha-1 everywhere -> anchors 10 Tg per cell.
- input fertilizer (2020 -> 2025): `regone.1` 100 -> 100, `regone.2`
  100 -> 150, `regtwo.1` 10 -> 10 (regional delta: `regone` +50, `regtwo` 0).
- projection puts the +50 on `regone.1` (only cell with area), `regone.2`'s
  10 Tg anchor hits zero area -> Inf -> 0 -> lost.

| quantity | current (measured) | expected after fix |
| --- | --- | --- |
| rates 2025 `regone.1` / `regone.2` / `regtwo.1` | 685.71 / 0 / **114.29** | 700 / 0 / **100** |
| region tonnages `regone` / `regtwo` | 68.57 / 11.43 | 70 / 10 |
| failed checks | none | none |

Derivation: regional expected `regone` = 20 + 50 = 70 over actual 60 ->
factor 7/6 -> 60 Tg/100 Mha * 7/6 = 700 kg ha-1; `regtwo` factor 1. Current
global factor = 80/70 scales `regtwo.1` to 114.29 -> the test is RED today on
the `regtwo.1 == 100` and `regone.1 == 700` expectations. Rates stay well below
the 1200 cap, so `toolCapFertilizerKgPerHa` is a no-op here.
Also assert `expect_false(any(grepl("Check failed", conditions$conditions)))`
and the unchanged target reproduction for 2010/2020.

**Fixture 2 — global fallback with warn check** (cells `regone.1`,
`regtwo.1`; `regone.1` -> 0 Mha in 2025; flat input; anchors 10 + 10 Tg):

- expected: `regone` 10 but actual 0 -> unachievable -> fallback: global
  factor 20/10 -> rates `regone.1` 0, `regtwo.1` **200** (measured: identical
  rates today, but WITHOUT any failed check).
- assert `expect_message(..., "Check failed.*cropland left to hold it")`
  (RED today: no such message) and additionally a failed
  `"Regional fertilizer totals are parallel"` check, while the GLOBAL check
  stays green.

Note: fixture names must keep the region prefix scheme `regone.N` /
`regtwo.N` — magclass splits dim-1 names at the dot, and test 1's
`"reg.one"`/`"reg.two"` are both region `reg` (that is exactly why existing
test 1 is unaffected by this change: regional == global there).

## Steps (TDD)

1. **Red tests** — append the two `test_that` blocks to
   `tests/testthat/test-calcNonlandHarmonized.R`, copying the fixtures from
   `/p/tmp/R_tmp/opencode/fert-regional/draft.R` and replacing the prints with
   the expectations from the table above. Confirm red with
   `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-calcNonlandHarmonized.R")'`
   (fixture 1 fails on 114.29 != 100 and 685.71 != 700; fixture 2 fails on the
   missing message).
2. **Green** — implement the code above in `R/calcNonlandHarmonized.R`
   (rescale block + regional check + roxygen). Confirm the two new tests pass
   AND all pre-existing tests in the file stay green.
3. **Regression** — `make format`, then `make test`; all green.
4. **Production verify (recommended)** — cache auto-invalidates (the cached
   `calc*` body changed). Rerun `/p/tmp/R_tmp/opencode/fert-compare/runfert.R`
   and compare against `/p/tmp/R_tmp/opencode/fert-compare/fert_proj.rds`
   (previous run) with a small measure script in
   `/p/tmp/R_tmp/opencode/fert-regional/`. Expect: global parallelism check
   green, `toolCheckFertilizer` green, regional totals essentially unchanged
   (only relocations across regions from the old global factor move back),
   the new warn check ideally never fires (if it fires, that is news —
   investigate, do not silently accept).
5. **Summary** — append a short paragraph to `executive-summary.md`: what the
   global factor did, what the regional rescale changes, measured relocation
   before-vs-after from step 4.

## Settled decisions / assumptions

- Scope: `absoluteChanges` only; fade/offset chains untouched.
- Fallback semantics (decision 1a + amend): global factor for the whole table
  of one crop x year whenever ANY region x crop x year has expected > 0 but
  actual == 0; reported via `toolExpectTrue(..., falseStatus = "warn")`, NOT a
  plain `warning()` (maintainer amendment).
- Strengthened per-region parallelism check included (decision 2a); in the
  fallback case it reports a failed check by design.
- Per-region clamping of negative expected replaces the global clamp;
  `expectedGlobalTg <- dimSums(expectedRegionTg, 1)` keeps the global check
  exact in both branches.
- `factor == 1` where actual == 0 and expected == 0 (regions without the crop
  at all).
- Cap (`toolCapFertilizerKgPerHa`) and `toolFertilizerKgPerHa` unchanged; cap
  runs after the rescale and already preserves region sums, so the regional
  property survives both steps (except in the cap's own global fallback case —
  unchanged pre-existing behavior).
- madrat cache: the change lives inside the cached `calcNonlandHarmonized`
  body itself, so the cache key updates automatically; `toolExpectTrue` /
  `toolExpectLessDiff` come from `mstools` (`Depends`) — no cache-key gotcha.
- Existing tests verified safe: test 1 single region `reg` (regional ==
  global), tests 3/4 have no compensation loss (factor 1 everywhere), the two
  `toolCapFertilizerKgPerHa` tests do not go through this code path.
- Never touch `man/`; do not run `make install`/`docs`/`build`.
- `if`/`for`/`while` always with `{` and newline.

## Environment note

`git` currently fails inside this worktree ("not a git repository:
/p/tmp/pascalfu/mrdownscale-absharm-full/.git/worktrees/fert-land-projection"
— stale worktree registration). If commits are requested, run
`git worktree repair` from the main repo
`/p/tmp/pascalfu/mrdownscale-absharm-full` first.
