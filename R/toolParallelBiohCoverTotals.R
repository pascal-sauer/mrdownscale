#' toolParallelBiohCoverTotals
#'
#' For the absoluteChanges harmonizer, keep the bioh cover totals (forest:
#' bioh.primf/secmf/secyf/pltns, other: bioh.primn/secnf) parallel to the input
#' (constant offset) after the harmonization year. The cover targets are
#' computed by toolParallelCoverTargets and allocated onto the recalculated
#' bioh by toolAllocateRegionTotals, split onto the cover's bioh items
#' proportional to their recalculated values.
#'
#' Processing steps, per year after the harmonization year:
#'
#' - Regional cover target: anchor from the target data in the harmonization
#'   year plus the input change relative to that year.
#' - Floor regional targets at zero and rescale the rest onto the raw global
#'   target, so regions with a loss target receive nothing while the global
#'   cover total stays input parallel.
#' - Floor a negative global cover target at zero and charge the clipped loss
#'   to the other cover, keeping the global total over all covers parallel; the
#'   regionally rescaled targets are the balanced targets, the floored global
#'   total the balanced global target.
#' - Allocate each cover's balanced regional target onto the recalculated bioh
#'   of its clusters (regional weight: summed recalculated bioh of the cover's
#'   items) and split it onto the cover's items proportional to their
#'   recalculated values, so clusters without harvest area host nothing.
#' - Check the allocation for regional and global cover parallelism to the
#'   balanced targets (toolExpectLessDiff, 10^-5).
#'
#' Output properties: In the ideal case, the bioh of
#' each cover stays parallel to the input after the harmonization year,
#' regionally and globally, and so does the sum over both covers. When that
#' is impossible, parallelism is broken as follows:
#'
#' 1. A regional target that turns negative is floored to zero and the
#'    remaining regions are rescaled, keeping the global cover total parallel;
#'    regional parallelism breaks.
#' 2. A cover whose global target turns negative is floored to zero and the
#'    loss charged to the other cover: that cover's global total is no longer
#'    parallel, the sum over both covers still is.
#' 3. Stranded targets of regions without recalculated bioh fill the remaining
#'    room of all clusters (toolAllocateRegionTotals), so receivers stay within
#'    their recalculated bioh weight; if the global room falls short, the
#'    unplaceable remainder is not allocated and breaks regional and global
#'    cover parallelism, flagged by the tool's own checks as notes (no error).
#' 4. Clusters without harvest area have zero recalculated bioh and host
#'    nothing.
#'
#' Always true: the allocated bioh of a cluster stays within its recalculated
#' bioh, except through the uncapped own share of an oversubscribed region,
#' and nothing lands on clusters without recalculated bioh.
#'
#' @param xInput recategorized input nonland data, containing the harmonization
#'   year and the years after it
#' @param xTarget extrapolated target nonland data, providing the anchor in the
#'   harmonization year, its years up to the harmonization year are returned
#'   unchanged
#' @param biohCalculated recalculated bioh (harmonized kg C per Mha times
#'   harmonized harvest area), providing the allocation weight, only its years
#'   after the harmonization year are used
#' @param harmonizationYear integer year the input changes are measured from
#' @return bioh for all years: the target data up to the harmonization year,
#'   then the allocated cover targets for the years after it; the tool checks
#'   its own allocation for regional and global parallelism of every cover to
#'   the balanced targets
#'
#' @author Pascal Sauer
toolParallelBiohCoverTotals <- function(xInput, xTarget, biohCalculated, harmonizationYear) {
  yearsAfterHy <- getYears(biohCalculated, as.integer = TRUE)
  yearsAfterHy <- yearsAfterHy[yearsAfterHy > harmonizationYear]
  covers <- list(forest = c("bioh.primf", "bioh.secmf", "bioh.secyf", "bioh.pltns"),
                 other = c("bioh.primn", "bioh.secnf"))
  inputCover <- list()
  targetCover <- list()
  coverItems <- list()
  for (cover in names(covers)) {
    items <- intersect(covers[[cover]], getItems(biohCalculated, 3))
    inputCover[[cover]] <- dimSums(xInput[, c(harmonizationYear, yearsAfterHy), items], c(1.2, 3))
    targetCover[[cover]] <- dimSums(xTarget[, harmonizationYear, items], c(1.2, 3))
    coverItems[[cover]] <- items
  }
  # the cover targets are floored, scaled and balanced such that the global
  # cover totals stay input parallel
  coverTotals <- toolParallelCoverTargets(inputCover, targetCover, harmonizationYear)
  targetCoverBioh <- coverTotals$coverTargets
  coverGlobalBalanced <- coverTotals$coverGlobalBalanced

  # allocate every cover's regional target onto the recalculated bioh
  # proportional to its weight (item split by the recalculated values), so
  # clusters without harvest area host nothing and the regional totals hold
  allocatedCoverBioh <- list()
  for (cover in names(targetCoverBioh)) {
    items <- coverItems[[cover]]
    coverBiohWeight <- dimSums(biohCalculated[, yearsAfterHy, items], 3)
    msgCover <- paste0("Regions with target ", cover, " bioh have harvest area left to hold it")
    itemShare <- biohCalculated[, yearsAfterHy, items] / coverBiohWeight
    itemShare[is.na(itemShare)] <- 0
    allocatedCoverBioh[[cover]] <-
      toolAllocateRegionTotals(targetCoverBioh[[cover]], coverBiohWeight, msgCover) * itemShare
  }

  # absoluteChanges must stay parallel to the input (constant offset) instead of matching it
  allocatedCoverTotals <- lapply(allocatedCoverBioh, function(allocated) dimSums(allocated, c(1.2, 3)))
  for (cover in names(allocatedCoverTotals)) {
    toolExpectLessDiff(allocatedCoverTotals[[cover]], targetCoverBioh[[cover]], 10^-5,
                       paste0("Regional ", cover, " bioh totals are parallel to input ",
                              "after the harmonization year"))
  }
  toolExpectLessDiff(dimSums(Reduce(`+`, allocatedCoverTotals), 1), coverGlobalBalanced, 10^-5,
                     "Global bioh totals are parallel to input after the harmonization year")

  biohAdapted <- do.call(mbind, allocatedCoverBioh)
  histYears <- getYears(xTarget, as.integer = TRUE)
  histYears <- histYears[histYears <= harmonizationYear]
  return(mbind(xTarget[, histYears, getItems(biohAdapted, 3)], biohAdapted))
}
