# Plan: project harmonized fertilizer tonnages onto the harmonized land (absoluteChanges)

## Context

Production run (`input = "magpie"`, `target = "luh3"`, `harmonizationPeriod =
c(2020, 2020)`, `harmonization = "absoluteChanges"`, full resolution) fails
`toolCheckFertilizer` ("Fertilizer application is <= 1200 kg ha-1 yr-1", max
1.3e10). Root cause: fertilizer tonnages are built from MAgPIE's land logic
but divided by the harmonized (LUH3-based) land (see below). This branch
(`fert-land-projection`, off `absharm-full` at b5775e1) fixes that directly.

NOTE: ignore the `fert-area-floor` branch. It is an abandoned, superseded
experiment (eps floor zeroing rates of cells with cropland < 10^-2 Mha; it
reduced the residue max to 22211 but left 62 cells / 28 cell x crop pairs
above 1200). Do NOT merge, cherry-pick, or reuse its code/tests — this branch
is self-contained and its projection + cap design fully replaces the floor
(the projection gives residue cells no delta by construction, and the cap
clamps and redistributes whatever anchor tonnage remains). Only its production
*measurements* are cited below as diagnostic context.

### Root cause (diagnosed, artifacts below)

Fertilizer is harmonized in **tonnage space**: per cell and crop group
`Tg = targetTg(hp) + (inputTg(t) - inputTg(hp))`
(`toolHarmonizeAbsoluteChanges`), then divided by **harmonized** cropland area
(`toolFertilizerKgPerHa`). The harmonized land is only ever used as a post-hoc
divisor — it never shapes the tonnes:

- Input tonnages are split into crop groups by **MAgPIE input land** shares
  (`toolRecategorizeFertilizer` in `calcNonlandInputRecategorized`), year by year.
- The target anchor uses LUH3 areas; harmonized land post-hp follows LUH3 at hp
  plus MAgPIE area *deltas*.

Step-0 diagnosis of the 28 pairs (`diag0.R`/`diag0.log`/`fert_diag0.rds` in
`/p/tmp/R_tmp/opencode/fert-compare/`, reconstruction matches pipeline tonnage
to 2e-16):

- **No cell-projection mismatch**: cell-total LUH3/input cropland ratios are
  0.76-1.39 for all pairs; datasets agree on cell totals.
- Disagreement is in the **crop-group composition inside cells** (group-area
  ratios 0.22-23, e.g. IND.52 c4ann: MAgPIE 4.95 vs LUH3 1.10 Mha at hp), and
  it is **dynamic**: 24/28 pairs have ~zero input fertilizer at hp (rate <= 30
  kg ha-1) and their tonnes are post-2020 MAgPIE intensification (IND.52 c4ann:
  0.011 -> 0.287 Tg by 2040), while the harmonized group **area collapses**
  along MAgPIE's own negative area trend (IND.52 c4ann: 1.10 -> 0.013 Mha;
  NEU.108 c3ann: 1.64 -> 0.031 Mha). Tons grow where the land truth says the
  crop vanishes. The remaining 4/28 are already-intensive cells (input rates
  up to 784 kg ha-1, LUH3 target rates up to 1015) pushed over 1200 by growth.
- A static hp-area reprojection cannot fix this: the mismatch grows over time.

## Design: full projection + residual cap

Decision (maintainer): treat the harmonized land as **truth** for fertilizer
tonnages in the `absoluteChanges` branch of `calcNonlandHarmonized`, in Tg
space, for years after the harmonization year only:

```
regionDelta(g, t) = regionSum(inputTg(g, t)) - regionSum(inputTg(g, hp))
Tg_h(cell, g, t)  = targetTg(hp, cell, g)                              # anchor, per-cell, unchanged
                  + regionDelta(g, t) * area_h(cell, g, t) / regionSum(area_h(., g, t))
```

The input delta (intensification) is redistributed within each region x group x
year proportional to the harmonized crop-group area of that year, so extra
tonnes land only where the truth says the crop is. Region sums are preserved;
the existing global parallelism rescale then runs unchanged. After the rescale
a single-pass **residual cap** handles anchor residues (e.g. NEU.108 c3ann
~1680 kg ha-1 with zero delta, and LUH3 area residues down to 1e-12 Mha): per
region x group x year,

```
excess_i = max(0, tg_i - 1200 * area_i / 10^3)      # donors land exactly at the cap
spare_i  = max(0, 1200 * area_i / 10^3 - tg_i)      # receivers stay below it
tg_i'    = tg_i - excess_i + (sum(excess) / sum(spare)) * spare_i
```

Headroom-weighting makes one vectorized pass sufficient (receivers scale their
own headroom, so no overshoot); **no iteration**. If `sum(excess) > sum(spare)`
in a region (infeasible there, in theory), fall back to one pass over the
**global** pool — allowed because the parallelism check is global; per step 0
this never triggers (max post-rescale region average ~884 <= 1200, headroom
huge: all 28 pairs' region averages <= 479 kg ha-1).

### Settled decisions / assumptions

- Scope: `absoluteChanges` only; fade/offset and all other chains untouched.
- Regional projection scope (preserves MAgPIE's between-region intensification
  structure; global scope would erase it).
- Weights are harmonized areas **of year t**, not hp (newly expanded cropland
  is itself part of the input growth and must receive its share).
- Negative `regionDelta` allowed; per-cell clamping at 0 loses a tiny amount
  which the existing global rescale absorbs (as today).
- The per-cell target anchor stays as-is (preserves LUH3 spatial pattern and
  pre-hp exact reproduction); the residual cap, not anchor projection, handles
  collapsed-area cells.
- No area floor is introduced: it would be redundant — residue cells get
  ~zero delta under the projection, and the cap clamps plus redistributes any
  remaining anchor tonnage (strictly subsuming the abandoned floor mechanism).
  The existing `Inf -> 0` / `is.nan -> 0` handlers for exactly-zero areas in
  `toolFertilizerKgPerHa` / `calcNonlandHarmonized` stay untouched.
- Cap threshold = 1200 (same constant as `toolCheckFertilizer`).

## Steps

1. **Red tests** — extend `tests/testthat/test-calcNonlandHarmonized.R`:
   - same-region fixture where one cell's harmonized crop-group area collapses
     after hp while its input tonnages grow: after harmonization its delta is
     zero and its rate <= 1200, the region's other cells absorbed the tonnes;
   - residue fixture: a cell with 10^-5 Mha harmonized cropland carrying target
     + input delta tonnes ends at rate <= 1200 with its excess redistributed
     to region partners (replaces the abandoned floor's `== 0` expectation);
   - region x group x year sums preserved by the projection;
   - global parallelism expectation still exact; pre-hp years still reproduce
     the target; non-fixture cells untouched.
   Confirm red.
2. **Green** — in the `absoluteChanges` block of `R/calcNonlandHarmonized.R`:
   replace the per-cell input delta with the region-projected delta above
   (helper code inlined in the cached body or a `tool*` function — madrat
   cache-key rule, AGENTS.md); add the single-pass residual cap as a new
   `tool*` function (auto cache-tracked), applied after the rescale; update
   the roxygen description. No changes to `toolFertilizerKgPerHa` or `man/`.
3. **Regression** — `make format`, then `make test`; all green.
4. **Production verify** — cache auto-invalidated; rerun
   `/p/tmp/R_tmp/opencode/fert-compare/runfert.R` and an exceedance-detail
   script (reuse `measure3.R` pattern). Expect: `toolCheckFertilizer` passes
   (0 cells > 1200 — including the former residue cells that reached 1.3e10
   without any floor), parallelism check green, pre-2020 reproduces target,
   rescale factors in the previous range. Plausibility: per-region group
   averages essentially unchanged vs before (projection preserves region sums;
   expect differences only from clamping, fractions of a Tg). Report global
   max rate and total relocated tonnes before vs after.
5. **Executive summary** — write `executive-summary.md` in this worktree root:
   short and concise — what was done, why, how, and how it improved the
   situation, with a before-vs-after table (cells > 1200: 192 -> 0, global max
   1.3e10 -> <= 1200, check: fails -> passes; measurement artifacts stay in
   /p/tmp/R_tmp/opencode/fert-compare/).

## Constraints

- Ignore the `fert-area-floor` branch entirely (abandoned; no merges,
  cherry-picks, or code reuse — see note in Context).
- TDD: red test first for every behavior (AGENTS.md).
- `tool*` helpers must be real `tool*` functions or inlined in the cached body
  (madrat cache-key gotcha).
- Pre-hp output must reproduce the target exactly (existing check).
- Global parallelism invariant: sums per crop group x year after hp stay
  parallel to input (`toolExpectLessDiff` 10^-5).
- Never touch `man/`; do not run `make install`/`docs`/`build`.
- If/for/while always with `{` and newline.
