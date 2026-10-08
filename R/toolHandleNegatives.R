#' toolHandleNegatives
#'
#' Negative values are set to zero and the remaining values are scaled down so
#' that the category sums of every cell and year match the target area, which
#' defaults to the category sums of x before zeroing.
#'
#' Categories listed in lastScaled are only scaled down after all other
#' categories have been scaled down to zero, so they are protected against
#' compensation cuts unless the other categories cannot absorb the excess area
#' alone. With the default lastScaled = character(0) all categories are scaled
#' down proportionally together.
#'
#' Intended to be applied to a group of land cover categories (e.g. forest
#' plus other land), so that negative values in one category are compensated
#' by the others without changing the total group area.
#'
#' @param x magclass object, may contain negative values
#' @param targetArea optional target area per cell and year, recycled if it has
#'   fewer years than x, must not be negative
#' @param lastScaled character vector of category names that are only scaled
#'   down after all other categories have been scaled down to zero, names not
#'   present in x are ignored
#' @return x without negative values, category sums matching the target area
#'
#' @author Pascal Sauer
toolHandleNegatives <- function(x, targetArea = NULL, lastScaled = character(0)) {
  tolerance <- 10^-5
  if (is.null(targetArea)) {
    targetArea <- dimSums(x, 3)
  }
  stopifnot(all(targetArea >= 0))
  x[x < 0] <- 0
  currentArea <- dimSums(x, 3)
  fact <- targetArea / (currentArea + (currentArea == 0))
  stopifnot(all(fact <= 1 + tolerance))
  fact[fact > 1] <- 1
  lastItems <- intersect(lastScaled, getItems(x, 3))
  if (length(setdiff(getItems(x, 3), lastItems)) == 0) {
    lastItems <- character(0)
  }
  if (length(lastItems) == 0) {
    x <- fact * x
  } else {
    excess <- pmax(currentArea - targetArea, 0)
    otherItems <- setdiff(getItems(x, 3), lastItems)
    otherArea <- dimSums(x[, , otherItems], 3)
    lastArea <- currentArea - otherArea
    otherFact <- (otherArea - pmin(excess, otherArea)) / (otherArea + (otherArea == 0))
    lastFact <- (lastArea - pmax(excess - otherArea, 0)) / (lastArea + (lastArea == 0))
    x[, , otherItems] <- otherFact * x[, , otherItems]
    x[, , lastItems] <- lastFact * x[, , lastItems]
  }
  stopifnot(all(abs(dimSums(x, 3) - targetArea) < tolerance))
  return(x)
}
