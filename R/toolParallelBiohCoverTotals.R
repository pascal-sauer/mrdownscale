#' toolParallelBiohCoverTotals
#'
#' For absoluteChanges, keep the regional forest and other land cover totals of
#' bioh parallel to the input (constant offset) after the harmonization year:
#' each cover's regional target is floored at zero and scaled onto the global
#' target, so a region with a negative change receives nothing while the rest
#' keep the global cover total parallel, and is then allocated onto the
#' recalculated bioh directly (toolAllocateRegionTotals), proportional to
#' harmonized harvest area times harvest rate, so clusters without harvest area
#' host nothing and a region with a bioh target but no harvest area left to
#' hold it hands its share to the global pool, where the receiving regions
#' exceed their own target, reported only by the warn-only checks. If the
#' global target of a cover itself turns negative, it is floored to zero and
#' the clipped loss is charged to the other covers by uniformly scaling their
#' targets, keeping the global total over all covers parallel.
#'
#' @param xInput recategorized input nonland data, containing the harmonization
#'   year and the years after it
#' @param xTarget extrapolated target nonland data, its years up to the
#'   harmonization year are returned unchanged
#' @param biohCalculated recalculated bioh (harmonized kg C per Mha times
#'   harmonized harvest area), only its years after the harmonization year are
#'   used
#' @param harmonizationYear integer year the input changes are measured from
#' @return bioh for all years: the target data up to the harmonization year,
#'   then the allocated cover targets for the years after it
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
