#' toolCapFertilizerKgPerHa
#'
#' Enforce a plausibility threshold on fertilizer rates by redistributing the
#' tonnage of every cluster above it over the remaining headroom of the
#' clusters of the same region, crop group and year, or of the global pool for
#' those crop group and year variables in which a region cannot compensate for its own
#' excess. Receivers only receive a share proportional to their own headroom,
#' so a single pass suffices and no one overshoots. Variables in which not even
#' the global pool has enough headroom are left unchanged, so that the tonnage
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
  # convert the threshold from kg ha-1 yr-1 to Tg yr-1 per cluster
  capTg <- threshold * cropMha * 10^-3
  cappedTg <- toolRedistributeOverCap(fertilizerTg, capTg,
                                      warnMessage = "Regions compensate for their own fertilizer above the threshold")
  cappedKgPerHa <- toolFertilizerKgPerHa(cappedTg, landMha)
  # rounding in the tonnage round trip can push values that landed exactly on
  # the cap slightly above it, the lost tonnes are far below all check
  # tolerances, infeasible variables stay unchanged for toolCheckFertilizer to report
  feasible <- dimSums(pmax(fertilizerTg - capTg, 0), 1) <= dimSums(pmax(capTg - fertilizerTg, 0), 1)
  roundedOver <- (cappedKgPerHa > threshold) * feasible
  cappedKgPerHa[roundedOver == 1] <- threshold
  return(cappedKgPerHa)
}
