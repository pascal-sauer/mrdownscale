#' toolParallelBiohCoverTotals
#'
#' For absoluteChanges, keep the regional forest and other land cover totals of
#' bioh parallel to the input (constant offset) after the harmonization year:
#' each cover's regional expectation is floored at zero and scaled onto the
#' global expectation, so a region expecting a loss receives nothing while the
#' rest keep the global cover total parallel, and is then allocated onto the
#' recalculated bioh directly (toolAllocateRegionTotals), proportional to
#' harmonized harvest area times harvest rate, so clusters without harvest area
#' host nothing and a region expecting bioh without any harvest area left to
#' hold it hands its share to the global pool, where the receiving regions
#' exceed their own expectation, reported only by the warn-only checks. If the
#' global expectation of a cover itself turns negative, it is floored to zero
#' and the clipped loss is charged to the other covers by uniformly scaling
#' their targets, keeping the global total over all covers parallel.
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
#'   then the allocated cover expectations for the years after it
#'
#' @author Pascal Sauer
toolParallelBiohCoverTotals <- function(xInput, xTarget, biohCalculated, harmonizationYear) {
  yearsAfterHy <- getYears(biohCalculated, as.integer = TRUE)
  yearsAfterHy <- yearsAfterHy[yearsAfterHy > harmonizationYear]
  covers <- list(forest = c("bioh.primf", "bioh.secmf", "bioh.secyf", "bioh.pltns"),
                 other = c("bioh.primn", "bioh.secnf"))
  expectedCoverBioh <- list()
  coverGlobal <- list()
  coverItems <- list()
  for (cover in names(covers)) {
    items <- intersect(covers[[cover]], getItems(biohCalculated, 3))
    if (length(items) == 0) {
      next
    }
    inputCoverBioh <- dimSums(xInput[, c(harmonizationYear, yearsAfterHy), items], c(1.2, 3))
    expected <- setYears(dimSums(xTarget[, harmonizationYear, items], c(1.2, 3)), NULL) +
      (inputCoverBioh[, yearsAfterHy, ] - setYears(inputCoverBioh[, harmonizationYear, ], NULL))
    # global first allocation: floor the regional cover expectations at zero
    # and scale the rest onto the raw global expectation, keeping the global
    # cover total parallel even where regional expectations had to be floored
    globalExpected <- dimSums(expected, 1)
    positiveTotal <- dimSums(pmax(expected, 0), 1)
    scale <- pmax(globalExpected, 0) / positiveTotal
    scale[is.na(scale)] <- 0
    expectedCoverBioh[[cover]] <- pmax(expected, 0) * scale
    coverGlobal[[cover]] <- globalExpected
    coverItems[[cover]] <- items
  }

  # the negative global cover expectation floored to zero is charged to the
  # other covers, so the global total over all covers stays input parallel
  globalTotal <- Reduce(`+`, coverGlobal)
  positiveGlobalTotal <- Reduce(`+`, lapply(coverGlobal, function(globalExpected) pmax(globalExpected, 0)))
  balanceFactor <- pmax(globalTotal, 0) / positiveGlobalTotal
  balanceFactor[is.na(balanceFactor)] <- 1
  coverGlobalBalanced <- pmax(globalTotal, 0)
  for (cover in names(expectedCoverBioh)) {
    expectedCoverBioh[[cover]] <- expectedCoverBioh[[cover]] * balanceFactor
  }

  # allocate every cover's regional expectation directly onto the
  # recalculated bioh, proportional to its weight (split onto the items by
  # their recalculated values, which together is proportional to them), so
  # clusters without harvest area host nothing and the regional cover totals
  # hold by construction
  allocatedCoverBioh <- list()
  for (cover in names(expectedCoverBioh)) {
    items <- coverItems[[cover]]
    coverBiohWeight <- dimSums(biohCalculated[, yearsAfterHy, items], 3)
    msgCover <- paste0("Regions with expected ", cover, " bioh have harvest area left to hold it")
    itemShare <- biohCalculated[, yearsAfterHy, items] / coverBiohWeight
    itemShare[is.na(itemShare)] <- 0
    allocatedCoverBioh[[cover]] <-
      toolAllocateRegionTotals(expectedCoverBioh[[cover]], coverBiohWeight, msgCover) * itemShare
  }

  # absoluteChanges must stay parallel to the input (constant offset) instead of matching it
  allocatedCoverTotals <- lapply(allocatedCoverBioh, function(allocated) dimSums(allocated, c(1.2, 3)))
  for (cover in names(allocatedCoverTotals)) {
    toolExpectLessDiff(allocatedCoverTotals[[cover]], expectedCoverBioh[[cover]], 10^-5,
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
