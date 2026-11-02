#' toolParallelBiohCoverTotals
#'
#' For the absoluteChanges harmonizer, keep the bioh cover totals (forest:
#' bioh.primf/secmf/secyf/pltns, other: bioh.primn/secnf) parallel to the input
#' (constant offset) after the harmonization year.
#'
#' toolParallelCoverTargets computes the regional cover targets (absolute
#' amounts per region and year). toolAllocateRegionTotals distributes each
#' regional target over the clusters of the region, proportional to their
#' recalculated bioh. Within a cluster, the amount is split over the cover's
#' bioh items, proportional to their recalculated bioh. Clusters with zero
#' recalculated bioh get nothing.
#'
#' Processing steps, per year after the harmonization year:
#'
#' - Regional cover target: the value from the target data in the
#'   harmonization year plus the change in the input since that year.
#' - Set negative regional targets to zero. This raises their sum. To
#'   compensate, scale the regional targets > 0 down, so the global cover
#'   total equals the raw global target. Regions set to zero get nothing.
#' - Set a negative global cover target to zero. The clipped loss is
#'   compensated for by the other cover: its regional targets are scaled by
#'   one common factor, so the global total over both covers stays input
#'   parallel. The scaled targets are the balanced targets, the resulting
#'   global total the balanced global target.
#' - Distribute each cover's balanced regional target over the clusters of the
#'   region, proportional to their recalculated bioh (regional weight: summed
#'   recalculated bioh of the cover's items), and split the amount of each
#'   cluster over the cover's items, proportional to their recalculated values.
#' - Check that the regional and global cover totals match the balanced
#'   targets (toolExpectLessDiff, 10^-5).
#'
#' Output properties: Ideally, the bioh of each cover stays parallel to the
#' input after the harmonization year, regionally and globally, and so does
#' the sum over both covers. When this is not possible, parallelism breaks as
#' follows:
#'
#' 1. A negative regional target is set to zero and the targets > 0 of the
#'    other regions are scaled down. The global cover total stays parallel,
#'    the parallelism of that region is lost.
#' 2. A cover whose global target is negative is set to zero, and the clipped
#'    loss is compensated for by the other cover. The global total of the
#'    cover set to zero is no longer parallel, the sum over both covers is.
#' 3. If a region has zero recalculated bioh, but target > 0, that bioh is
#'    moved to clusters in other regions (toolAllocateRegionTotals). If there
#'    is not enough capacity left, regional and global cover parallelism
#'    breaks, flagged by the tool's own checks as notes (no error).
#'
#' Always true: a cluster never gets more than its recalculated bioh. The one
#' exception is a region whose target exceeds its total recalculated bioh:
#' there, every cluster gets its proportional share of the target, even if it
#' exceeds its own recalculated bioh. Clusters with zero recalculated bioh
#' never get anything.
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
  # cover targets: negative values set to zero, scaled and balanced so the
  # global cover totals stay input parallel
  coverTotals <- toolParallelCoverTargets(inputCover, targetCover, harmonizationYear)
  targetCoverBioh <- coverTotals$coverTargets
  coverGlobalBalanced <- coverTotals$coverGlobalBalanced

  # distribute every cover's regional target over the clusters of the region
  # (weight: recalculated bioh, item split by the recalculated values), so
  # clusters without harvest area get nothing and the regional totals hold
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
  globalAllocated <- 0
  for (cover in names(allocatedCoverTotals)) {
    toolExpectLessDiff(allocatedCoverTotals[[cover]], targetCoverBioh[[cover]], 10^-5,
                       paste0("Regional ", cover, " bioh totals are parallel to input ",
                              "after the harmonization year"))
    globalAllocated <- globalAllocated + allocatedCoverTotals[[cover]]
  }
  toolExpectLessDiff(dimSums(globalAllocated, 1), coverGlobalBalanced, 10^-5,
                     "Global bioh totals are parallel to input after the harmonization year")

  biohAdapted <- do.call(mbind, allocatedCoverBioh)
  histYears <- getYears(xTarget, as.integer = TRUE)
  histYears <- histYears[histYears <= harmonizationYear]
  return(mbind(xTarget[, histYears, getItems(biohAdapted, 3)], biohAdapted))
}
