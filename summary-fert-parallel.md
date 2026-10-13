# Executive summary: global fertilizer parallel to input under absoluteChanges

## What

The `absoluteChanges` fertilizer lines in `plotHarmonizedNonLand()` (global totals,
Tg yr-1) now run **parallel to the unharmonized input lines from the harmonization
year on**, per fertilizer crop category. Years up to the harmonization year still
return the target data unchanged.

## Why

The raw absolute-changes formula `target(hy) + input(t) - input(hy)` is parallel in
global sums by construction, but two cell-local corrections broke that afterwards:
zero clamping of negative harmonized values (clamped cells stick at 0 while the
input declines, pushing the global sum above parallel) and the Inf->0 replacement
where harmonized cropland vanished but raw fertilizer Tg > 0 (fertilizer vanishes,
pulling the global sum below parallel). Both deviations grew over time.

## How

`calcNonlandHarmonized()` (R/calcNonlandHarmonized.R) now rescales the fertilizer
kg ha-1 yr-1 rates per crop category and year after the Inf->0 handling, so the
implied global totals (rate x harmonized cropland, via `toolFertilizerTg`) equal
`max(0, targetGlobal(hy) + (inputGlobal(t) - inputGlobal(hy)))`. The scaling factor
is `expected / actual` (1 where `actual == 0`); since area is unchanged, scaling
the rate scales Tg by the same factor. This absorbs both clamping excess and
Inf->0 losses, redistributing fertilizer lost in cells without cropland to cells
that still have some. The helper computing the expected totals is defined inside
`calcNonlandHarmonized` (madrat cache-key rule). A `toolExpectLessDiff` check now
asserts the parallel global totals for years after the harmonization year,
replacing reliance on the plot; the input-match check for other harmonizations is
unchanged. `toolHarmonizeAbsoluteChanges` is untouched.

Verified with new unit tests (synthetic magclass data with mocked `calcOutput`,
including the Inf->0 redistribution case: a lost 140 Tg rescales the remaining
cell's rate by 290/150) and on real MAgPIE/LUH3 data via the cached pipeline:
max |delta(harmonized) - delta(input)| across all crop categories and years > 2020
is 5.7e-14 Tg. The extreme per-cell rates in near-zero-cropland cells (~10^10
kg ha-1) are pre-existing and unaffected by the rescale.

## Reproduce

```sh
Rscript -e 'suppressMessages(pkgload::load_all(".")); plotHarmonizedNonLand(harmonizations = c("absoluteChanges", "inputUnharmonized"))'
```

Writes `plot-harmonized-nonland-fertilizer.png` (plus bioh / harvest / wood harvest
area groups) to the working directory; the fertilizer panel shows the
absoluteChanges lines running parallel to the dashed input lines from the 2020
vline onwards. Offline guarantee: `make test` (179 tests green).
