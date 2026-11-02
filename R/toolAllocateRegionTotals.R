#' toolAllocateRegionTotals
#'
#' Distribute target amounts (e.g. wood harvest area) per region and variable
#' over the clusters (weight: their capacity, e.g. forest land), so the totals
#' of the regions with capacity above zero reach their target exactly. A
#' cluster can never physically receive more than its capacity. Each cluster
#' first receives its own share: the region's target distributed over the
#' region's clusters proportional to their capacity. Shares are not capped at
#' the cluster's capacity. A region whose target exceeds its region's total
#' capacity gives its clusters more than their capacity, on purpose, and the
#' caller has to account for the excess. Targets > 0 of regions without any
#' capacity are pooled globally and spread over the remaining room (capacity
#' minus own share) of all clusters, so borrowed amounts never push a cluster
#' past its capacity. If the global room is too small for the whole stranded
#' pool, the remainder is not allocated and has to be accounted for by the
#' caller.
#'
#' This is a helper function for harmonization, not a general tool to
#' disaggregate region totals over clusters.
#'
#' @param targetRegionAmount target amount total per region and variable, e.g.
#'   wood harvest area
#' @param clusterCapacity capacity per cluster, doubling as the distribution
#'   weight, e.g. forest land
#' @param warnMessage message of the warn check whether regions with target
#'   amounts > 0 have receiving capacity > 0
#' @return amounts per cluster, proportional to their capacity
#'
#' @author Pascal Sauer
toolAllocateRegionTotals <- function(targetRegionAmount, clusterCapacity,
                                     warnMessage = "Regions with target totals > 0 have receiving capacity") {
  # a capacity of zero can never receive anything, setting to 1 only avoids NaN
  regionCapacity <- dimSums(clusterCapacity, 1.2)
  hasCapacity <- regionCapacity > 0
  regionCapacity[regionCapacity == 0] <- 1
  toolExpectTrue(all(hasCapacity | targetRegionAmount == 0), warnMessage, falseStatus = "warn")
  # the own target share of every cluster within its region, uncapped, so a
  # region beyond its total capacity overshoots its clusters on purpose
  ownShare <- targetRegionAmount * hasCapacity / regionCapacity * clusterCapacity
  # the pooled stranded targets of regions without capacity fill the room
  # (capacity minus own share) left on all clusters, unallocated if it falls short
  room <- pmax(clusterCapacity - ownShare, 0)
  strandedFill <- dimSums(targetRegionAmount * (1 - hasCapacity), 1) / dimSums(room, 1)
  strandedFill[is.na(strandedFill)] <- 0
  strandedFill[strandedFill > 1] <- 1
  return(ownShare + room * strandedFill)
}
