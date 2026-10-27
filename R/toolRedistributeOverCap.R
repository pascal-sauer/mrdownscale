#' toolRedistributeOverCap
#'
#' Enforce a cap on per-cell amounts by redistributing the amount of every cell
#' above it onto the remaining headroom of the cells of the same region and
#' table, or of the global pool for those tables in which a region cannot absorb
#' its own excess. Receivers only receive a share proportional to their own
#' headroom, so a single pass suffices and no one overshoots. Tables in which
#' not even the global pool has enough headroom are left unchanged, so that the
#' amount totals stay intact and the caller can report the exceedance.
#'
#' @param xAmount amounts per cell
#' @param capAmount cap per cell, in the same unit as xAmount
#' @param message message of the warn check whether regions absorb their own excess
#' @return amounts per cell, capped, with the totals of the feasible tables preserved
#'
#' @author Pascal Sauer
toolRedistributeOverCap <- function(xAmount, capAmount, message = "Regions absorb their own excess") {
  excess <- pmax(xAmount - capAmount, 0)
  spare <- pmax(capAmount - xAmount, 0)

  # infeasible tables stay unchanged, capping would lose amounts, the caller reports them
  globalExcess <- dimSums(excess, 1)
  globalSpare <- dimSums(spare, 1)
  feasible <- globalExcess <= globalSpare
  excess <- excess * feasible
  spare <- spare * feasible
  globalExcess <- globalExcess * feasible

  # per table, the global pool only where a region cannot absorb its own excess
  regionExcess <- dimSums(excess, 1.2)
  regionSpare <- dimSums(spare, 1.2)
  useGlobal <- dimSums(regionExcess > regionSpare, 1) > 0
  toolExpectTrue(!any(useGlobal), message, falseStatus = "warn")
  # where spare is zero, excess is zero too or the other pool applies, so 1 only avoids NaN
  regionSpare[regionSpare == 0] <- 1
  globalSpare[globalSpare == 0] <- 1
  factor <- regionExcess / regionSpare * (1 - useGlobal) + globalExcess / globalSpare * useGlobal

  return(xAmount - excess + factor * spare)
}
