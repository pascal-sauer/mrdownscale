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
#' - Regional targets are set to zero where negative and the rest is rescaled
#'   onto the raw global target, so regions with a loss target receive nothing
#'   while the global cover total stays input parallel.
#' - Set a negative global cover target to zero and compensate the clipped loss
#'   with the other cover, keeping the global total over all covers parallel.
#'   The regional targets rescaled as a result of this compensation are the
#'   balanced targets, their global sums the balanced global target.
#' - Keep the primary harvest land-derived and distribute the remaining target
#'   over the clusters (weight: their harvestable secondary land), so the
#'   physical yearly cap for primary harvest holds by construction.
#' - Move harvest above the own cover's cap into the headroom of the other
#'   cover, region first and then globally. Only what fits nowhere becomes
#'   excess.
#'
#' Output properties: In the ideal case, the harvest area of
#' each cover stays parallel to the input after the harmonization year,
#' regionally and globally, and so does the sum over both covers. When that
#' is impossible, parallelism is broken as follows:
#'
#' 1. A regional target that turns negative is set to zero and the
#'    remaining regions are rescaled, keeping the global cover total parallel;
#'    regional parallelism breaks.
#' 2. A cover whose global target turns negative is set to zero and the loss
#'    is compensated for by the other cover: that cover's global total is no
#'    longer parallel, the sum over both covers still is.
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
#'   does not fit even in the other cover (`excessSecdHarvest`); the tool
#'   checks its own allocation for regional and global parallelism of every
#'   cover to the input
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
  coverTarget <- list()
  coverGlobal <- list()
  coverItems <- list()
  coverTotal <- list()
  coverCap <- list()
  coverLoss <- list()

  for (cover in names(covers)) {
    prim <- covers[[cover]]$prim
    secd <- covers[[cover]]$secd
    inputCover <- dimSums(inputRegion[, , c(prim, secd)], 3)
    targetCover <- dimSums(targetRegion[, , c(prim, secd)], 3)
    target <- setYears(targetCover[, hy, ], NULL) +
      (inputCover[, futureYears, ] - setYears(inputCover[, hy, ], NULL))
    # global first allocation: set negative regional targets to zero and
    # scale the rest onto the raw global target, keeping the global cover
    # total input parallel even where regional targets had to be set to zero
    globalTarget <- dimSums(target, 1)
    positiveTotal <- dimSums(pmax(target, 0), 1)
    scale <- pmax(globalTarget, 0) / positiveTotal
    scale[is.na(scale)] <- 0
    coverTarget[[cover]] <- pmax(target, 0) * scale
    coverGlobal[[cover]] <- globalTarget
  }

  # the negative global cover target set to zero is compensated for by the
  # other cover, so the global total over all covers stays input parallel
  globalTotal <- coverGlobal$forest + coverGlobal$other
  positiveGlobalTotal <- pmax(coverGlobal$forest, 0) + pmax(coverGlobal$other, 0)
  balanceFactor <- pmax(globalTotal, 0) / positiveGlobalTotal
  balanceFactor[is.na(balanceFactor)] <- 1
  coverGlobalBalanced <- pmax(globalTotal, 0)
  for (cover in names(covers)) {
    coverTarget[[cover]] <- coverTarget[[cover]] * balanceFactor
  }

  for (cover in names(covers)) {
    secd <- covers[[cover]]$secd
    requiredSecd <- coverTarget[[cover]] - dimSums(dimSums(primHarv[, , covers[[cover]]$prim], 1.2), 3)
    requiredSecd[requiredSecd < 0] <- 0

    # distribute each region's cover target over its clusters (weight: their
    # harvestable land, item split proportional to the item caps), the
    # clamped remainder above the cap flowing into the other cover's headroom
    msgHold <- paste0("Regions with target ", cover, " harvest area > 0 have harvestable land left to hold it")
    maxHarvItems <- maxHarv[, , secd]
    clusterCap <- dimSums(maxHarvItems, 3)
    itemShare <- maxHarvItems / clusterCap
    itemShare[is.na(itemShare)] <- 0
    clusterLevelAlloc <- toolAllocateRegionTotals(requiredSecd, clusterCap, msgHold)
    # stranded target that fits nowhere in the remaining cluster room is
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

  # harvest above the own cap is shifted into the headroom of the other cover,
  # region first and then globally, only what fits nowhere is reported as
  # excess
  for (cover in names(covers)) {
    other <- setdiff(names(covers), cover)
    for (margin in c(1.2, 1)) {
      excess <- pmax(coverTotal[[cover]] - coverCap[[cover]], 0)
      headroomItems <- pmax(maxHarv[, , covers[[other]]$secd] - coverItems[[other]], 0)
      headroom <- dimSums(headroomItems, 3)
      fitShare <- dimSums(headroom, margin) / dimSums(excess, margin)
      fitShare[is.na(fitShare)] <- 0
      fitShare[fitShare > 1] <- 1
      filled <- dimSums(excess, margin) / dimSums(headroom, margin)
      filled[is.na(filled)] <- 0
      filled[filled > 1] <- 1
      coverTotal[[cover]] <- coverTotal[[cover]] - excess * fitShare
      coverTotal[[other]] <- coverTotal[[other]] + headroom * filled
      # the clusters of the other cover receive the transferred harvest within
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
