#' toolCapFertilizerKgPerHa
#'
#' Enforce a plausibility threshold on fertilizer rates by redistributing the
#' tonnage of every cell above it onto the remaining headroom of the cells of
#' the same region, crop group and year, or of the global pool for those crop
#' group and year tables in which a region cannot absorb its own excess.
#' Receivers only receive a share proportional to their own headroom, so a
#' single pass suffices and no one overshoots. Tables in which not even the
#' global pool has enough headroom are left unchanged, so that the tonnage
#' totals stay intact and toolCheckFertilizer reports the exceedance.
#'
#' @param fertilizerKgPerHa fertilizer rates in kg ha-1 yr-1
#' @param landMha land in Mha, used to convert to and from tonnages
#' @param threshold maximum plausible rate in kg ha-1 yr-1
#' @return fertilizer rates in kg ha-1 yr-1, capped at the threshold
#'
#' @author Pascal Sauer
toolCapFertilizerKgPerHa <- function(fertilizerKgPerHa, landMha, threshold = 1200) {
  cropMha <- toolAggregateCropland(landMha[getItems(fertilizerKgPerHa, 1), getItems(fertilizerKgPerHa, 2), ],
                                   getItems(fertilizerKgPerHa, 3.2), keepOthers = FALSE)
  fertilizerTg <- toolFertilizerTg(fertilizerKgPerHa, landMha)
  # convert the threshold from kg ha-1 yr-1 to Tg yr-1 per cell
  capTg <- threshold * cropMha * 10^-3
  excessTg <- pmax(fertilizerTg - capTg, 0)
  spareTg <- pmax(capTg - fertilizerTg, 0)

  # infeasible tables stay unchanged, capping would lose tonnes, toolCheckFertilizer reports them
  globalExcessTg <- dimSums(excessTg, 1)
  globalSpareTg <- dimSums(spareTg, 1)
  feasible <- globalExcessTg <= globalSpareTg
  excessTg <- excessTg * feasible
  spareTg <- spareTg * feasible
  globalExcessTg <- globalExcessTg * feasible

  # per crop group and year, the global pool only where a region cannot absorb its own excess
  regionExcessTg <- dimSums(excessTg, 1.2)
  regionSpareTg <- dimSums(spareTg, 1.2)
  useGlobal <- dimSums(regionExcessTg > regionSpareTg, 1) > 0
  toolExpectTrue(!any(useGlobal), "Regions absorb their own fertilizer above the threshold", falseStatus = "warn")
  # where spare is zero, excess is zero too or the other pool applies, so 1 only avoids NaN
  regionSpareTg[regionSpareTg == 0] <- 1
  globalSpareTg[globalSpareTg == 0] <- 1
  factor <- regionExcessTg / regionSpareTg * (1 - useGlobal) + globalExcessTg / globalSpareTg * useGlobal

  cappedKgPerHa <- toolFertilizerKgPerHa(fertilizerTg - excessTg + factor * spareTg, landMha)
  # rounding in the tonnage round trip can push values that landed exactly on
  # the cap slightly above it, the lost tonnes are far below all check
  # tolerances
  roundedOver <- (cappedKgPerHa > threshold) * feasible
  cappedKgPerHa[roundedOver == 1] <- threshold
  return(cappedKgPerHa)
}
