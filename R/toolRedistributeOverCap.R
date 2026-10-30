#' toolRedistributeOverCap
#'
#' Enforce a cap on per-cluster amounts by redistributing the amount of every
#' cluster above it onto the remaining headroom of the clusters of the same
#' region and variable, or of the global pool for those variables in which a
#' region cannot absorb its own excess. Receivers only receive a share
#' proportional to their own headroom, so a single pass suffices and no one
#' overshoots. Variables in which not even the global pool has enough headroom
#' are left unchanged, so that the amount totals stay intact and the caller can
#' report the exceedance.
#'
#' @param xAmount amounts per cluster
#' @param capAmount cap per cluster, in the same unit as xAmount
#' @param warnMessage message of the warn check whether regions absorb their own excess
#' @return amounts per cluster, capped, with the totals of the feasible variables preserved
#'
#' @author Pascal Sauer
toolRedistributeOverCap <- function(xAmount, capAmount, warnMessage = "Regions absorb their own excess") {
  excess <- pmax(xAmount - capAmount, 0)
  spare <- pmax(capAmount - xAmount, 0)

  # infeasible variables stay unchanged, capping would lose amounts, the caller reports them
  globalExcess <- dimSums(excess, 1)
  globalSpare <- dimSums(spare, 1)
  feasible <- globalExcess <= globalSpare
  excess <- excess * feasible
  spare <- spare * feasible
  globalExcess <- globalExcess * feasible

  # per variable, the global pool only where a region cannot absorb its own excess
  regionExcess <- dimSums(excess, 1.2)
  regionSpare <- dimSums(spare, 1.2)
  useGlobal <- dimSums(regionExcess > regionSpare, 1) > 0
  toolExpectTrue(!any(useGlobal), warnMessage, falseStatus = "warn")
  # where spare is zero, excess is zero too or the other pool applies, so 1 only avoids NaN
  regionSpare[regionSpare == 0] <- 1
  globalSpare[globalSpare == 0] <- 1
  factor <- regionExcess / regionSpare * (1 - useGlobal) + globalExcess / globalSpare * useGlobal

  return(xAmount - excess + factor * spare)
}
