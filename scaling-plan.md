# Plan: protect primf/primn and urban in the negative-value rescaling of `toolHarmonizeAbsoluteChanges`

Follow-up to PR review comment on `R/toolHarmonizeAbsoluteChanges.R` (thread
`PRRT_kwDOLdnbrM6hd0wi`, "is 0 <= factor <= 1 always true here?").

## Problem

In the clamp block (`if (any(changed < 0))`), after negative values are set to 0 the
non-prim categories are scaled by
`factor <- (targetArea - primSum) / dimSums(changed[, , nonPrim], dim = 3)`
to keep the cell total at `targetArea` (total target area in the harmonization year).

Two issues:

1. `urban` is currently scaled first-class together with secdf/secdn/etc., although it
   should be protected.
2. `factor < 0` (or `-Inf`) occurs exactly when `primSum > targetArea`, i.e. prim
   overshoots the total because prim is never rescaled while non-forest categories went
   negative and were clipped to 0 (e.g. per cell `primf = 90, primn = 20, urban = -5,
   other = -5`, `targetArea = 100` after clipping: prim = 110 > 100). The current code
   only warns and clamps the factor to 0, which silently emits too much total area
   (since the interim commit that added the final consistency `stopifnot`, this case
   makes the tool **error** instead).

## Required behavior (per review discussion)

Reduce the cell total from `T = primSum + urbanSum + otherSum >= targetArea` back to
`targetArea` with this priority:

1. Scale all remaining categories (`other` = everything except primf/primn/urban, i.e.
   pltns, secdf, secdn, cropland, ...) as far as needed, down to 0 at most.
2. Only if prim + urban still exceed the target area afterwards, scale primf/primn
   proportionally (e.g. factor = 10/11 for prim = 110, targetArea = 100).
3. Only if urban alone exceeds the target area, scale urban as last resort.

Each step's factor is then guaranteed to be in [0, 1]. Note `toolReplaceExpansion`
afterwards still caps primf/primn at the previous year's level, so any remaining prim
gain shows up as secdf/secdn expansion anyway.

## Implementation

Replace the body of `if (any(changed < 0))` (keep `changed[changed < 0] <- 0`; the
`targetArea` definition was hoisted to the top of the function in the meantime):

```r
    prim <- intersect(c("primf", "primn"), getItems(changed, dim = 3))
    urban <- intersect("urban", getItems(changed, dim = 3))
    other <- setdiff(getItems(changed, dim = 3), c(prim, urban))
    sumOf <- function(items) if (length(items) > 0) dimSums(changed[, , items], dim = 3) else 0
    primSum <- sumOf(prim)
    urbanSum <- sumOf(urban)
    otherSum <- sumOf(other)
    # 1. scale remaining categories (never prim/urban) as far as needed and possible
    factor <- (targetArea - primSum - urbanSum) / (otherSum + (otherSum == 0))
    factor[!is.finite(factor) | factor < 0] <- 0
    if (length(other) > 0) changed[, , other] <- changed[, , other] * pmin(factor, 1)
    # 2. if prim + urban still exceed the target area, scale prim
    if (any(primSum + urbanSum > targetArea + 10^-5)) {
      toolStatusMessage("warn", paste0("prim + urban area exceed total area after correcting ",
                                       "negative values, scaling prim categories down"),
                        level = level)
    }
    factor <- (targetArea - urbanSum) / (primSum + (primSum == 0))
    factor[!is.finite(factor) | factor < 0] <- 0
    if (length(prim) > 0) changed[, , prim] <- changed[, , prim] * pmin(factor, 1)
    # 3. urban alone exceeding the total area: scale urban as last resort
    if (any(urbanSum > targetArea + 10^-5)) {
      toolStatusMessage("warn", paste0("urban area exceeds total area after correcting negative ",
                                       "values, scaling urban down"),
                        level = level)
    }
    factor <- targetArea / (urbanSum + (urbanSum == 0))
    factor[!is.finite(factor) | factor < 0] <- 0
    if (length(urban) > 0) changed[, , urban] <- changed[, , urban] * pmin(factor, 1)
```

Mechanics:

- `pmin(factor, 1)`: `factor < 1` exactly when the still-protected tiers already exceed
  `targetArea`; unused tiers get `factor = 1` and stay untouched.
- `+(sum == 0)` guards the division by 0 (same idiom as `safeAvailable` in the forest
  deduction block); where the tier sum is 0 all its values are 0, so multiplying by any
  finite factor is a no-op.
- Upper bound of each factor follows from the precondition `dimSums(changed, dim = 3) >=
  targetArea - 10^-5` (clipping only raises cell sums, and before clipping they equal
  `targetArea`, enforced by the `stopifnot` further up); `pmin(..., 1)` absorbs the
  float slack of the 10^-5 tolerance.
- The old `factor < 0` warning inside the block is replaced by the two new tier
  warnings.
- Final consistency check at the end of the function (`total == targetArea`) now passes
  in all cases instead of erroring on prim overshoot.

## Also update

- Roxygen description: replace "are set to 0 and the non-prim categories are then
  scaled, so that the total area stays constant." with wording that states the
  protection order: other categories first, primf/primn only if prim + urban exceed the
  total area, urban as last resort.
- `man/toolHarmonizeAbsoluteChanges.Rd` via `devtools::document()`.

## Tests (`tests/testthat/test-toolHarmonizeAbsoluteChanges.R`)

1. Existing test "clamps and scales non-prim categories if other land is insufficient"
   must be updated: with urban protected, `reg.five` 2025 becomes
   `(primf 40, primn 0, secdf 0, secdn 0, urban 5, other 55)` instead of
   `(40, 0, 0, 0, 60/13, 720/13)`. Rename the test accordingly (e.g. "scales remaining
   categories while protecting prim and urban").
2. New test, prim scaling: target 2010/2020 `(40, 4, 10, 1, 5, 40)`; input 2020
   `(20, 4, 20, 1, 5, 50)`, 2025 `(90, 10, 0, 0, 0, 0)` (all totals 100) ->
   `changed = (110, 10, -10, 0, 0, -10)`; forest deduction zeroes secdf/primn/secdn,
   clipping zeroes other, prim = 110 > 100 -> prim scaled by 10/11 -> primf 100;
   `toolReplaceExpansion` moves the 40 -> 100 primf expansion into secdf, so expect
   2025 = `(40, 0, 60, 0, 0, 0)`, total 100, expect prim warning.
3. New test, urban last resort: target 2010/2020 `(0, 0, 10, 5, 80, 5)`; input 2020
   `(0, 0, 10, 5, 10, 75)`, 2025 `(0, 0, 0, 0, 100, 0)` (totals 100) ->
   `changed = (0, 0, 0, 0, 170, -70)`; clipping zeroes other, urban 170 > 100 -> urban
   scaled by 100/170 -> 100; expect 2025 urban 100, total 100, expect urban warning.
4. All tests keep asserting `dimSums(out, dim = 3) == 100` and `all(out >= 0)`.

## PR reply draft (thread PRRT_kwDOLdnbrM6hd0wi)

Fixed range explanation plus your overshoot point, with the refinement from the
follow-up discussion: factor <= 1 always held for the scaling tier (cell sums equal
`targetArea` before clipping via the new total-area `stopifnot`; clipping only increases
them), while `factor < 0` occurred exactly when the protected categories overshot
`targetArea`. Instead of clamping to 0 (silently emitting too much area), the
correction now has an explicit priority: first scale all non-prim/non-urban categories
down (to 0 max), scale primf/primn proportionally (e.g. 10/11 for prim = 110 vs total
100) only if prim + urban still exceed the total, and scale urban only if urban alone
exceeds the total. Every factor is therefore guaranteed to be in [0, 1], and the total
area is conserved in all cases.
