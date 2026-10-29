#' toolAllocateRegionTotals
#'
#' Allocate expected amounts per region and table directly onto the cells,
#' proportional to their weight, so the totals of the regions with hosting
#' capacity (weight above zero) reach their expectation exactly. Expectations
#' of regions without any capacity are pooled globally and distributed over
#' all cells proportional to their weight as well.
#'
#' @param expectedRegionAmount expected amount total per region and table
#' @param cellWeight weight per cell, e.g. the available hosting capacity,
#'   cells with weight zero never receive anything
#' @param message message of the warn check whether regions with expected
#'   amounts have hosting capacity
#' @return amounts per cell, proportional to their weight, with the totals of
#'   the regions with hosting capacity matching their expectation plus their
#'   share of the globally pooled expectations of the regions without any
#'
#' @author Pascal Sauer
toolAllocateRegionTotals <- function(expectedRegionAmount, cellWeight,
                                     message = "Regions with expected totals have hosting capacity") {
  receiver <- dimSums(cellWeight, 1.2) > 0
  # a weight of zero can never receive anything, the 1 only avoids NaN
  regionWeight <- dimSums(cellWeight, 1.2)
  regionWeight[regionWeight == 0] <- 1
  globalWeight <- dimSums(cellWeight, 1)
  globalWeight[globalWeight == 0] <- 1
  toolExpectTrue(all(receiver | expectedRegionAmount == 0), message, falseStatus = "warn")
  globalStranded <- dimSums(expectedRegionAmount * (1 - receiver), 1)
  return(((expectedRegionAmount * receiver / regionWeight) + globalStranded / globalWeight) * cellWeight)
}
