# Plan: global fertilizer parallel to input under absoluteChanges (Option A)

Worktree `.worktrees/fert-parallel`, branch `fert-parallel`, base `absharm-full` @ 8109a76.
Do NOT touch `man/`, do NOT run `make install/docs/build`. No commits unless asked.

## Goal

In `plotHarmonizedNonLand()` (global totals, Tg yr-1), the `absoluteChanges` fertilizer
lines must be **parallel to the unharmonized input lines from the harmonization year on**
(per fertilizer crop category, i.e. per plotted line). Level offset relative to input is
fine and intended; only the year-to-year deltas must match.

## Diagnosis (why it is not parallel today)

Fertilizer is harmonized in Tg yr-1 in `R/calcNonlandHarmonized.R:40-61` via
`toolHarmonizeAbsoluteChanges(constantTotal = FALSE)`:
`raw = targetTg(hy) + (inputTg(t) - inputTg(hy))`.
The global sum of this unclamped formula is exactly parallel (constant offset
`targetGlobalTg(hy) - inputGlobalTg(hy)`). Parallelism is broken afterwards by:

1. per-cell zero clamping (`R/toolHarmonizeAbsoluteChanges.R:116`): clamped cells stick
   at 0 while input declines -> global harmonized drifts **above** parallel;
2. Inf->0 loss in `R/calcNonlandHarmonized.R:54-61`: where harmonized cropland was
   clamped to 0 Mha but raw fertilizer Tg > 0, kg/ha rate is Inf and set to 0 -> Tg
   vanishes -> global falls **below** parallel;
3. both are cell-local, so the deviation grows over time in both directions.

## Chosen approach: Option A (keep clamping + global rescale)

Keep `toolHarmonizeAbsoluteChanges` untouched. In `calcNonlandHarmonized`, inside the
`harmonization == "absoluteChanges"` branch and AFTER the Inf->0 handling, rescale the
fertilizer kg/ha rates per crop category x year so the implied global Tg totals
(rate x harmonized cropland, via `toolFertilizerTg`) equal:

- years <= hy: unchanged (target data, leave as is);
- years > hy: `max(0, targetGlobalTg(hy) + (inputGlobalTg(t) - inputGlobalTg(hy)))`
  per crop category.

Notes:

- `factor = expected / actual`, guard `actual == 0` (then factor 1);
- scaling kg/ha rates by the factor scales Tg by the same factor (same area);
- this absorbs both clamping excess and Inf->0 losses -> exact global parallelism;
- put the helper INSIDE the `calcNonlandHarmonized` body (madrat cache-key rule, see
  AGENTS.md); no new tool file needed;
- add a `toolExpectLessDiff` check in `calcNonlandHarmonized` asserting the parallel
  global totals for years > hy (replaces reliance on the plot);
- update the stale comment at `R/calcNonlandHarmonized.R:112-115` (deviation from input
  remains as a constant offset, deltas are now parallel) and the roxygen description;
  the `harmonization != "absoluteChanges"` input-match check stays as is.

## Steps (TDD)

1. Write failing tests in `tests/testthat/test-calcNonlandHarmonized.R` (synthetic
   magclass data, `local_mocked_bindings(calcOutput = ...)`, offline):
   - assert per-year global harmonized fertilizer Tg (rate x harmonized cropland) ==
     `targetGlobalTg(hy) + inputGlobalDelta` for all years > hy;
   - Inf->0 case (reg.two 2025 losing 140 Tg) must now be redistributed to the cell
     with remaining cropland, preserving the global parallel total (reg.one rate
     scales by 290/150 -> 2900 kg/ha in that fixture);
   - confirm red: `Rscript -e 'pkgload::load_all("."); testthat::test_file("tests/testthat/test-calcNonlandHarmonized.R")'`
2. Implement in `R/calcNonlandHarmonized.R` as sketched above; confirm green.
   Beware: rescaled rates near/above the 1200 kg/ha `toolCheckFertilizer` threshold
   will emit messages -> adjust first test to `expect_message` or tweak fixture values
   so unchanged paths stay warning-free.
3. `make test` (all green), then `make format`.
4. Try `make plot-harmonized-nonland` in the worktree to verify visually; without
   madrat cache/mainpath (this sandbox) it will likely fail -> run it on the cluster
   instead. Expected: fertilizer PNG shows absoluteChanges lines running parallel to
   the dashed input lines from the 2020 vline onwards.
5. Final step: write a concise `executive-summary.md` in the worktree folder:
   - what was done, why, how (few sentences each, no filler);
   - copy-paste-ready Rscript call that runs part of the pipeline so the changes are
     visible in output and relevant plots are created, using existing plotting only
     (no custom plotting), e.g.:

     ```sh
     Rscript -e 'suppressMessages(pkgload::load_all(".")); plotHarmonizedNonLand(harmonizations = c("absoluteChanges", "inputUnharmonized"))'
     ```

     which writes `plot-harmonized-nonland-fertilizer.png` (plus bioh / harvest
     groups) to the working directory.

## Assumptions

- Parallel per fertilizer crop category `fertilizer.<crop>` (the plotted lines); there
  is no nutrient dimension (MAgPIE gdx fertilizer is already aggregated over N/P/K).
- Years <= harmonization year keep returning target data unchanged.
- Negative expected global totals floor at 0 (crop line zeroed for that year).
- Real-data plot regeneration is the maintainer's cluster run; unit tests are the
  offline guarantee.
