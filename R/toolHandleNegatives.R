#' toolHandleNegatives
#'
#' Negative values are set to zero and the remaining values are scaled down so
#' that the category sums of every cell and year match the target area, which
#' defaults to the category sums of x before zeroing.
#'
#' Intended to be applied to a group of land cover categories (e.g. forest
#' plus other land), so that negative values in one category are compensated
#' by the others without changing the total group area.
#'
#' @param x magclass object, may contain negative values
#' @param targetArea optional target area per cell and year, recycled if it has
#'   fewer years than x
#' @param allowNegativeTarget if TRUE, categories of cells with a negative
#'   target area are set to zero instead of raising an error
#' @return x without negative values, category sums matching the target area
#'
#' @author Pascal Sauer
toolHandleNegatives <- function(x, targetArea = NULL, allowNegativeTarget = FALSE) {
  tolerance <- 10^-5
  if (is.null(targetArea)) {
    targetArea <- dimSums(x, 3)
  }
  if (allowNegativeTarget) {
    targetArea[targetArea < 0] <- 0
  } else {
    stopifnot(all(targetArea >= 0))
  }
  x[x < 0] <- 0
  currentArea <- dimSums(x, 3)
  fact <- targetArea / (currentArea + (currentArea == 0))
  stopifnot(all(-tolerance <= fact & fact <= 1 + tolerance))
  fact[fact > 1] <- 1
  x <- fact * x
  stopifnot(all(abs(dimSums(x, 3) - targetArea) < tolerance))
  return(x)
}
