#' toolRescaleRegionTotals
#'
#' Rescale per-cell amounts so that each region's amount total of every table
#' reaches its expected amount, distributing amounts lost where cells could no
#' longer hold them to the other cells of the same region, falling back to the
#' global factor exactly in those tables where a region expects amounts but has
#' no capacity left to hold any.
#'
#' @param xAmount amounts per cell
#' @param expectedRegionAmount expected amount total per region and table
#' @param message message of the warn check whether regions with expected
#'   amounts have hosting capacity
#' @return amounts per cell, rescaled so the region totals of the tables without
#'   global fallback match their expectation
#'
#' @author Pascal Sauer
toolRescaleRegionTotals <- function(xAmount, expectedRegionAmount,
                                    message = "Regions with expected totals have hosting capacity") {
  actualRegion <- dimSums(xAmount, 1.2)
  factor <- expectedRegionAmount / actualRegion
  factor[actualRegion == 0] <- 1
  # unholdable amounts: expected > 0 where the region holds nothing, zero
  # elsewhere, summed up per table
  unholdable <- dimSums(expectedRegionAmount * (actualRegion == 0), 1)
  if (!toolExpectTrue(sum(unholdable[]) == 0, message, falseStatus = "warn")) {
    # a region that holds nothing can never hold its amounts, so exactly the
    # affected tables fall back to the global factor
    actualGlobal <- dimSums(actualRegion, 1)
    factorGlobal <- dimSums(expectedRegionAmount, 1) / actualGlobal
    factorGlobal[actualGlobal == 0] <- 1
    fallback <- unholdable > 0
    factor <- factor * (fallback != TRUE) + factorGlobal * fallback
  }
  return(xAmount * factor)
}
