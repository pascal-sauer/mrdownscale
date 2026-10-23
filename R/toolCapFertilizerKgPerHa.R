#' toolCapFertilizerKgPerHa
#'
#' Enforce a plausibility threshold on fertilizer rates by redistributing the
#' tonnage of every cell above it onto the remaining headroom of the cells of
#' the same region, crop group and year, or of the global pool if a region
#' cannot absorb its own excess. Receivers only receive a share proportional to
#' their own headroom, so a single pass suffices and no one overshoots. If not
#' even the global pool has enough headroom, the input is returned unchanged so
#' that the tonnage totals stay intact and toolCheckFertilizer reports the
#' exceedance.
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
  regionExcessTg <- dimSums(excessTg, 1.2)
  regionSpareTg <- dimSums(spareTg, 1.2)
  if (any(regionExcessTg[] > regionSpareTg[])) {
    # at least one region cannot absorb its own excess, run the pass over the
    # global pool, whose totals are what the parallelism check constrains
    globalExcessTg <- dimSums(excessTg, 1)
    globalSpareTg <- dimSums(spareTg, 1)
    if (any(globalExcessTg[] > globalSpareTg[])) {
      # infeasible, capping would lose tonnes, let toolCheckFertilizer report it
      return(fertilizerKgPerHa)
    }
    factor <- globalExcessTg / globalSpareTg
  } else {
    # only needed where spare and excess are both zero, which yields factor 0
    regionSpareTg[regionSpareTg == 0] <- 1
    factor <- regionExcessTg / regionSpareTg
  }
  cappedKgPerHa <- toolFertilizerKgPerHa(fertilizerTg - excessTg + factor * spareTg, landMha)
  # rounding in the tonnage round trip can push values that landed exactly on
  # the cap slightly above it, the lost tonnes are far below all check
  # tolerances
  cappedKgPerHa[cappedKgPerHa > threshold] <- threshold
  return(cappedKgPerHa)
}
