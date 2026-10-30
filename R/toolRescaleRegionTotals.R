#' toolRescaleRegionTotals
#'
#' Rescale per-cluster amounts so that each region's amount total of every
#' variable reaches its target amount, distributing amounts lost where clusters
#' could no longer hold them to the other clusters of the same region, falling
#' back to the global factor exactly in those variables where a region expects
#' amounts but has no capacity left to hold any.
#'
#' @param xAmount amounts per cluster
#' @param targetRegionAmount target amount total per region and variable
#' @param warnMessage message of the warn check whether regions with target
#'   amounts have hosting capacity
#' @return amounts per cluster, rescaled so the region totals of the variables
#'   without global fallback match their target
#'
#' @author Pascal Sauer
toolRescaleRegionTotals <- function(xAmount, targetRegionAmount,
                                    warnMessage = "Regions with target totals have hosting capacity") {
  actualRegion <- dimSums(xAmount, 1.2)
  factor <- targetRegionAmount / actualRegion
  factor[actualRegion == 0] <- 1
  # unholdable amounts: target > 0 where the region holds nothing, zero
  # elsewhere, summed up per variable
  unholdable <- dimSums(targetRegionAmount * (actualRegion == 0), 1)
  if (!toolExpectTrue(sum(unholdable[]) == 0, warnMessage, falseStatus = "warn")) {
    # a region that holds nothing can never hold its amounts, so exactly the
    # affected variables fall back to the global factor
    actualGlobal <- dimSums(actualRegion, 1)
    factorGlobal <- dimSums(targetRegionAmount, 1) / actualGlobal
    factorGlobal[actualGlobal == 0] <- 1
    fallback <- unholdable > 0
    factor <- factor * (fallback != TRUE) + factorGlobal * fallback
  }
  return(xAmount * factor)
}
