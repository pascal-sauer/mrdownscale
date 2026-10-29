#' toolAllocateRegionTotals
#'
#' Allocate target amounts (e.g. wood harvest area) per region and variable
#' directly onto the clusters, proportional to their weight
#' (e.g. forest land), so the totals of the regions with hosting capacity
#' (weight above zero) reach their target exactly. Targets > 0 of regions
#' without any capacity are pooled globally and distributed over all clusters
#' proportional to their weight as well. For example, the target amount can be
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
  globalWeight <- dimSums(clusterWeight, 1)
  globalWeight[globalWeight == 0] <- 1
  toolExpectTrue(all(receiver | targetRegionAmount == 0), warnMessage, falseStatus = "warn")
  # the own target spread within the region plus the pooled stranded targets
  # of regions without capacity spread over all clusters
  regionalShare <- targetRegionAmount * receiver / regionWeight
  strandedShare <- dimSums(targetRegionAmount * (1 - receiver), 1) / globalWeight
  return((regionalShare + strandedShare) * clusterWeight)
}
