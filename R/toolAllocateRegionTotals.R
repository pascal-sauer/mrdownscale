#' toolAllocateRegionTotals
#'
#' Allocate target amounts (e.g. wood harvest area) per region and variable
#' directly onto the clusters, proportional to their weight
#' (e.g. forest land), so the totals of the regions with hosting capacity
#' (weight above zero) reach their target exactly. Own shares are not capped,
#' an oversubscribed region intentionally overshoots its clusters. Targets > 0
#' of regions without any capacity are pooled globally and spread over the
#' remaining room (weight minus own share) of all clusters, so no cluster ever
#' receives more than its weight on top of its own share. If the global room
#' cannot host the whole stranded pool, the remainder is not allocated and has
#' to be accounted for by the caller. For example, the target amount can be
#' the wood harvest area per region, the weight / capacity the forest land per
#' cluster which can host it.
#'
#' @param targetRegionAmount target amount total per region and variable, e.g.
#'   wood harvest area
#' @param clusterWeight weight/capacity per cluster, e.g. forest land
#' @param warnMessage message of the warn check whether regions with target
#'   amounts > 0 have hosting capacity/weight > 0
#' @return amounts per cluster, proportional to their weight
#'
#' @author Pascal Sauer
toolAllocateRegionTotals <- function(targetRegionAmount, clusterWeight,
                                     warnMessage = "Regions with target totals > 0 have hosting capacity") {
  # a weight of zero can never receive anything, setting to 1 only avoids NaN
  regionWeight <- dimSums(clusterWeight, 1.2)
  receiver <- regionWeight > 0
  regionWeight[regionWeight == 0] <- 1
  toolExpectTrue(all(receiver | targetRegionAmount == 0), warnMessage, falseStatus = "warn")
  # the own target share of every cluster within its region, uncapped so an
  # oversubscribed region overshoots on purpose
  ownShare <- targetRegionAmount * receiver / regionWeight * clusterWeight
  # the pooled stranded targets of regions without capacity fill the room
  # (weight minus own share) left on all clusters, unallocated if it falls short
  room <- pmax(clusterWeight - ownShare, 0)
  strandedFill <- dimSums(targetRegionAmount * (1 - receiver), 1) / dimSums(room, 1)
  strandedFill[is.na(strandedFill)] <- 0
  strandedFill[strandedFill > 1] <- 1
  return(ownShare + room * strandedFill)
}
