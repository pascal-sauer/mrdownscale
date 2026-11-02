#' toolParallelCoverTargets
#'
#' The cover target is the absolute cover per region and year: the value from
#' the target data in the harmonization year plus the change in the input
#' since that year.
#'
#' Per cover, the tool sets negative regional targets to zero. This raises
#' their sum. To compensate, it scales the targets > 0 of the same cover in
#' the other regions down, so the global total equals the raw (unmodified)
#' global target. Covers do not interact here.
#'
#' The tool also sets a negative global cover target to zero. It compensates
#' for the clipped loss by scaling the regional targets of all covers with one
#' common factor, so the global total over all covers is met (zero if it is
#' negative). Returns the balanced regional cover targets and the balanced
#' global cover total.
#'
#' @param inputCover named list of region-level per-year cover totals, the
#'   input including the harmonization year and the years after it
#' @param targetCover named list of region-level per-year cover totals,
#'   providing the anchor in the harmonization year
#' @param harmonizationYear integer year the input changes are measured from
#' @return list with the balanced regional cover targets per cover
#'   (`coverTargets`) and the balanced global cover total
#'   (`coverGlobalBalanced`)
#'
#' @author Pascal Sauer
toolParallelCoverTargets <- function(inputCover, targetCover, harmonizationYear) {
  yearsAfterHy <- getYears(inputCover[[1]], as.integer = TRUE)
  yearsAfterHy <- yearsAfterHy[yearsAfterHy > harmonizationYear]
  coverTarget <- list()
  globalTotal <- 0
  positiveGlobalTotal <- 0
  for (cover in names(inputCover)) {
    target <- setYears(targetCover[[cover]][, harmonizationYear, ], NULL) +
      (inputCover[[cover]][, yearsAfterHy, ] -
         setYears(inputCover[[cover]][, harmonizationYear, ], NULL))
    # set negative regional targets to zero and scale the rest down onto the
    # raw global target, keeping the global cover total input parallel
    globalTarget <- dimSums(target, 1)
    positiveTotal <- dimSums(pmax(target, 0), 1)
    scale <- pmax(globalTarget, 0) / positiveTotal
    scale[is.na(scale)] <- 0
    coverTarget[[cover]] <- pmax(target, 0) * scale
    globalTotal <- globalTotal + globalTarget
    positiveGlobalTotal <- positiveGlobalTotal + pmax(globalTarget, 0)
  }
  # compensate for a global cover target set to zero by scaling all regional
  # targets, so the global total over all covers stays input parallel
  balanceFactor <- pmax(globalTotal, 0) / positiveGlobalTotal
  balanceFactor[is.na(balanceFactor)] <- 1
  coverGlobalBalanced <- pmax(globalTotal, 0)
  coverTargets <- lapply(coverTarget, function(target) target * balanceFactor)
  return(list(coverTargets = coverTargets,
              coverGlobalBalanced = coverGlobalBalanced))
}
