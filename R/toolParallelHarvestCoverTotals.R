#' toolParallelHarvestCoverTotals
#'
#' For absoluteChanges, keep the global forest and other land cover totals of
#' the wood harvest area parallel to the input (constant offset) after the
#' harmonization year: the regional targets are floored at zero and scaled onto
#' the global target, so a region whose target is a loss receives nothing while
#' the remaining regions keep the global total parallel. If the global target of
#' a cover itself turns negative, it is floored to zero and the clipped loss is
#' charged to the other cover, keeping the global total over all covers
#' parallel. The prim harvest stays land-derived, the secd/pltns harvest of a
#' cover is allocated onto the clusters directly, proportional to their
#' harvestable secondary land, so a region with enough land meets its cover
#' target exactly and the physical cap holds by construction. Harvest that the
#' own cover cannot host fills the headroom of the other cover, region first
#' and then globally, only what remains unhostable anywhere is reported as
#' excess. Where the land-derived prim harvest alone outweighs a cover's
#' balanced target, the global totals may still exceed anchor plus delta, an
#' accepted physical residual reported only by the warn-only checks of the
#' caller.
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
#'   after the harmonization year (`secdHarv`), the excess harvest area that
#'   not even the other cover can host (`excessSecdHarvest`), the balanced
#'   regional cover targets (`coverTarget`), the balanced global cover total
#'   target (`coverGlobalBalanced`) and the cover structure (`covers`), so the
#'   caller can check its final output for parallelism
#'
#' @author Pascal Sauer
toolParallelHarvestCoverTotals <- function(xInput, xTarget, primHarv, landHarmonized, harmonizationYear) {
  hp <- harmonizationYear
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

  for (cover in names(covers)) {
    prim <- covers[[cover]]$prim
    secd <- covers[[cover]]$secd
    inputCover <- dimSums(inputRegion[, , c(prim, secd)], 3)
    targetCover <- dimSums(targetRegion[, , c(prim, secd)], 3)
    target <- setYears(targetCover[, hp, ], NULL) +
      (inputCover[, futureYears, ] - setYears(inputCover[, hp, ], NULL))
    # global first allocation: floor the regional targets at zero and
    # scale the rest onto the raw global target, keeping the global cover
    # total input parallel even where regional targets had to be floored
    globalTarget <- dimSums(target, 1)
    positiveTotal <- dimSums(pmax(target, 0), 1)
    scale <- pmax(globalTarget, 0) / positiveTotal
    scale[is.na(scale)] <- 0
    coverTarget[[cover]] <- pmax(target, 0) * scale
    coverGlobal[[cover]] <- globalTarget
  }

  # the negative global cover target floored to zero is charged to the
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

    # allocate the cover target of every region directly onto its clusters,
    # proportional to their harvestable land (split onto the items by their
    # caps, which together is proportional to the item caps themselves), so
    # the cap holds by construction and the clamped away remainder above the
    # cap flows on into the headroom of the other cover below
    msgHold <- paste0("Regions with target ", cover, " harvest area have harvestable land left to hold it")
    maxHarvItems <- maxHarv[, , secd]
    clusterCap <- dimSums(maxHarvItems, 3)
    itemShare <- maxHarvItems / clusterCap
    itemShare[is.na(itemShare)] <- 0
    allocatedSecd <- toolAllocateRegionTotals(requiredSecd, clusterCap, msgHold) * itemShare
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
    pmax(coverTotal$other - coverCap$other, 0)

  return(list(secdHarv = secdHarv,
              excessSecdHarvest = excessSecdHarvest,
              coverTarget = coverTarget,
              coverGlobalBalanced = coverGlobalBalanced,
              covers = covers))
}
