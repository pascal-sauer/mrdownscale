#' toolParallelCoverTargets
#'
#' Keep the global cover totals parallel to the input: floor every cover's
#' regional target at zero and rescale the remaining regions onto the raw
#' global target, then floor a negative global cover target at zero, charging
#' the clipped loss to the other covers by uniformly scaling their targets.
#' Returns the balanced regional cover targets and the balanced global cover
#' total.
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
  coverGlobal <- list()
  for (cover in names(inputCover)) {
    target <- setYears(targetCover[[cover]][, harmonizationYear, ], NULL) +
      (inputCover[[cover]][, yearsAfterHy, ] -
         setYears(inputCover[[cover]][, harmonizationYear, ], NULL))
    # floor the regional targets at zero and rescale the rest onto the raw
    # global target, keeping the global cover total input parallel
    globalTarget <- dimSums(target, 1)
    positiveTotal <- dimSums(pmax(target, 0), 1)
    scale <- pmax(globalTarget, 0) / positiveTotal
    scale[is.na(scale)] <- 0
    coverTarget[[cover]] <- pmax(target, 0) * scale
    coverGlobal[[cover]] <- globalTarget
  }
  # charge a negative global cover target floored to zero to the other covers,
  # so the global total over all covers stays input parallel
  globalTotal <- Reduce(`+`, coverGlobal)
  positiveGlobalTotal <- Reduce(`+`, lapply(coverGlobal, function(globalTarget) pmax(globalTarget, 0)))
  balanceFactor <- pmax(globalTotal, 0) / positiveGlobalTotal
  balanceFactor[is.na(balanceFactor)] <- 1
  coverGlobalBalanced <- pmax(globalTotal, 0)
  coverTargets <- lapply(coverTarget, function(target) target * balanceFactor)
  return(list(coverTargets = coverTargets,
              coverGlobalBalanced = coverGlobalBalanced))
}
