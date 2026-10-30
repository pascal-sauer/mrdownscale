#' toolParallelHarvestCoverTotals
#'
#' For the absoluteChanges harmonizer, keep the global wood harvest area totals
#' of the two covers (forest: primf/secdf/pltns, other: primn/secdn) parallel to
#' the input (constant offset) after the harmonization year.
#'
#' Processing steps, per year after the harmonization year:
#'
#' - Regional cover target: anchor from the target data in the harmonization
#'   year plus the input change relative to that year.
#' - Floor regional targets at zero and rescale the rest onto the raw global
#'   target, so regions with a loss target receive nothing while the global
#'   cover total stays input parallel.
#' - Floor a negative global cover target at zero and charge the clipped loss
#'   to the other cover, keeping the global total over all covers parallel; the
#'   regional targets rescaled by this charging are the balanced targets, their
#'   global sums the balanced global target.
#' - Keep the primary harvest land-derived and allocate the remaining target
#'   onto the clusters proportional to their harvestable secondary land, so the
#'   physical yearly cap holds by construction.
#' - Move harvest that the own cover cannot host into the headroom of the other
#'   cover, region first and then globally; only what no cover can host
#'   anywhere becomes excess.
#'
#' Output properties: In the ideal case, the harvest area of
#' each cover stays parallel to the input after the harmonization year,
#' regionally and globally, and so does the sum over both covers. When that
#' is impossible, parallelism is broken as follows:
#'
#' 1. A regional target that turns negative is floored to zero and the
#'    remaining regions are rescaled, keeping the global cover total parallel;
#'    regional parallelism breaks.
#' 2. A cover whose global target turns negative is floored to zero and the
#'    loss charged to the other cover: that cover's global total is no longer
#'    parallel, the sum over both covers still is.
#' 3. Where harvestable land is too little, harvest moves into the other
#'    cover's headroom, region first and then globally: first regional per
#'    cover parallelism is lost, then global per cover parallelism, the sum
#'    over both covers holds.
#' 4. The sum over both covers breaks where land-derived primary
#'    harvest alone outweighs a cover's balanced target, harvest exceeds
#'    anchor plus delta.
#' 5. Stranded cover targets of regions without any harvestable land for that
#'    cover fill the remaining room of all clusters (toolAllocateRegionTotals);
#'    what does not fit there is reported as excess directly, not shiftable.
#'
#' Always true: the allocated secondary and plantation harvest stays
#' within the harvestable land of every cluster.
#'
#' @param xInput recategorized wood harvest area input data, containing the
#'   harmonization year and the years after it
#' @param xTarget extrapolated wood harvest area target data, providing the
#'   anchor in the harmonization year
#' @param primHarv land-derived primary harvest area for the years after the
#'   harmonization year, with the items "primf" and "primn"
#' @param landHarmonized harmonized land data, providing the harvestable
#'   secondary land
#' @param harmonizationYear integer year the input changes are measured from
#' @return list with the secondary and plantation harvest area of the years
#'   after the harmonization year (`secdHarv`) and the excess harvest area that
#'   not even the other cover can host (`excessSecdHarvest`); the tool checks
#'   its own allocation for regional and global parallelism of every cover to
#'   the input
#'
#' @author Pascal Sauer
toolParallelHarvestCoverTotals <- function(xInput, xTarget, primHarv, landHarmonized, harmonizationYear) {
  hy <- harmonizationYear
  futureYears <- getYears(primHarv, as.integer = TRUE)
  inputRegion <- dimSums(toolAggregateWoodHarvest(xInput), 1.2)
  targetRegion <- dimSums(toolAggregateWoodHarvest(xTarget), 1.2)
  maxHarv <- toolMaxHarvestPerYear(landHarmonized, split = FALSE)[, futureYears, ]
  # the harvest covers, forest and other land, each with its prim source and
  # its secondary/plantation pool, target and allocation work per cover
  covers <- list(forest = list(prim = "primf", secd = c("secdf", "pltns")),
                 other = list(prim = "primn", secd = "secdn"))
  inputCover <- list()
  targetCover <- list()
  coverItems <- list()
  coverTotal <- list()
  coverCap <- list()
  coverLoss <- list()

  for (cover in names(covers)) {
    prim <- covers[[cover]]$prim
    secd <- covers[[cover]]$secd
    inputCover[[cover]] <- dimSums(inputRegion[, , c(prim, secd)], 3)
    targetCover[[cover]] <- dimSums(targetRegion[, , c(prim, secd)], 3)
  }
  # floor the regional cover targets at zero, rescale the rest onto the raw
  # global target and charge floored negative global cover targets to the
  # other cover, keeping the global total over all covers input parallel
  coverTotals <- toolParallelCoverTargets(inputCover, targetCover, hy)
  coverTarget <- coverTotals$coverTargets
  coverGlobalBalanced <- coverTotals$coverGlobalBalanced

  for (cover in names(covers)) {
    secd <- covers[[cover]]$secd
    requiredSecd <- coverTarget[[cover]] - dimSums(dimSums(primHarv[, , covers[[cover]]$prim], 1.2), 3)
    requiredSecd[requiredSecd < 0] <- 0

    # allocate each region's cover target onto its clusters proportional to
    # their harvestable land (item split proportional to the item caps), the
    # clamped remainder above the cap flowing into the other cover's headroom
    msgHold <- paste0("Regions with target ", cover, " harvest area > 0 have harvestable land left to hold it")
    maxHarvItems <- maxHarv[, , secd]
    clusterCap <- dimSums(maxHarvItems, 3)
    itemShare <- maxHarvItems / clusterCap
    itemShare[is.na(itemShare)] <- 0
    clusterLevelAlloc <- toolAllocateRegionTotals(requiredSecd, clusterCap, msgHold)
    # stranded target that no remaining cluster room could host anywhere is
    # pure loss, reported as excess proportional to the cluster caps
    # (uniformly where no cluster has any cap)
    unplaced <- pmax(dimSums(requiredSecd, 1) - dimSums(clusterLevelAlloc, 1), 0)
    lossShare <- clusterCap / dimSums(clusterCap, 1)
    lossShare[is.na(lossShare)] <- 1 / length(getCells(clusterCap))
    coverLoss[[cover]] <- lossShare * unplaced
    allocatedSecd <- clusterLevelAlloc * itemShare
    coverTotal[[cover]] <- dimSums(allocatedSecd, 3)
    coverCap[[cover]] <- clusterCap
    coverItems[[cover]] <- pmin(allocatedSecd, maxHarvItems)
  }

  # harvest of one cover that its own harvestable land cannot host is shifted
  # into the headroom of the other cover, region first and then globally,
  # only what remains unhostable anywhere is reported as excess
  for (cover in names(covers)) {
    other <- setdiff(names(covers), cover)
    for (margin in c(1.2, 1)) {
      excess <- pmax(coverTotal[[cover]] - coverCap[[cover]], 0)
      headroomItems <- pmax(maxHarv[, , covers[[other]]$secd] - coverItems[[other]], 0)
      headroom <- dimSums(headroomItems, 3)
      hostable <- dimSums(headroom, margin) / dimSums(excess, margin)
      hostable[is.na(hostable)] <- 0
      hostable[hostable > 1] <- 1
      filled <- dimSums(excess, margin) / dimSums(headroom, margin)
      filled[is.na(filled)] <- 0
      filled[filled > 1] <- 1
      coverTotal[[cover]] <- coverTotal[[cover]] - excess * hostable
      coverTotal[[other]] <- coverTotal[[other]] + headroom * filled
      # the transferred harvest lands on the clusters of the other cover within
      # their caps, proportional to their remaining headroom
      coverItems[[other]] <- coverItems[[other]] + headroomItems * filled
    }
  }

  secdHarv <- mbind(coverItems$forest, coverItems$other)
  excessSecdHarvest <- pmax(coverTotal$forest - coverCap$forest, 0) +
    pmax(coverTotal$other - coverCap$other, 0) +
    coverLoss$forest + coverLoss$other

  # self check of the final allocation: land-derived primary harvest plus the
  # allocated cover items stay parallel to the input, regionally and globally
  regionalCoverTotals <- list()
  for (cover in names(covers)) {
    regionalCoverTotals[[cover]] <- dimSums(dimSums(primHarv[, , covers[[cover]]$prim], 1.2), 3) +
      dimSums(dimSums(coverItems[[cover]], 1.2), 3)
    msg <- paste0("Regional ", cover, " wood harvest area is parallel to input ",
                  "after the harmonization year")
    toolExpectLessDiff(regionalCoverTotals[[cover]], coverTarget[[cover]], 10^-5, msg)
  }
  toolExpectLessDiff(dimSums(regionalCoverTotals$forest + regionalCoverTotals$other, 1),
                     coverGlobalBalanced, 10^-5,
                     "Global wood harvest area is parallel to input after the harmonization year")

  return(list(secdHarv = secdHarv,
              excessSecdHarvest = excessSecdHarvest))
}
