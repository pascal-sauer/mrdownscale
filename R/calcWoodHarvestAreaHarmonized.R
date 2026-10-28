#' calcWoodHarvestAreaHarmonized
#'
#' Harmonize wood harvest area based on harmonized land data.
#' First, wood harvest area is harmonized just like land data.
#' Then, to ensure consistency, primary harvest is converted to seondary
#' harvest (or vice versa) to match the primary land reduction from the
#' harmonized land data. If secondary harvest exceeds the available
#' secondary land, the excess is shifted to the other secondary land
#' (secdf to secdn and vice versa). Remaining excess harvest area is reported.
#' For absoluteChanges, the global forest and other land cover totals of the
#' harvest area stay parallel to the input (constant offset) after the
#' harmonization year: the regional expectations are floored at zero and scaled
#' onto the global expectation, so a region expecting a loss receives nothing
#' while the remaining regions keep the global total parallel. If the global
#' expectation of a cover itself turns negative, it is floored to zero and the
#' clipped loss is charged to the other cover, keeping the global total over all
#' covers parallel. The prim harvest
#' stays land-derived, the secd/pltns harvest of a cover absorbs the difference
#' within the region and the global pool of its own cover, falling back to the
#' headroom of the other cover, region first and then globally, and finally
#' reporting the excess where even that cannot host it.
#'
#' @inheritParams calcNonlandInput
#' @param target name of a target dataset, see \code{\link{calcLandTarget}}
#' for available target datasets
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param harmonization harmonization method, see \code{\link{toolGetHarmonizer}} for available methods
#' @return harmonized wood harvest area data
#'
#' @examples
#' \dontrun{
#' calcOutput("WoodHarvestAreaHarmonized", input = "magpie",
#'              target = "luh2mod", harmonizationPeriod = c(2015, 2050),
#'              harmonization = "fade")
#' }
#' @author Pascal Sauer
calcWoodHarvestAreaHarmonized <- function(input, target, harmonizationPeriod, harmonization) {
  landHarmonized <- calcOutput("LandHarmonized", input = input, target = target,
                               harmonizationPeriod = harmonizationPeriod,
                               harmonization = harmonization, aggregate = FALSE)

  # calculate prim harvest based on prim land reduction
  prim <- c("primf", "primn")
  primDiff <- setYears(landHarmonized[, -nyears(landHarmonized), prim],
                       getYears(landHarmonized)[-1]) - landHarmonized[, -1, prim]
  stopifnot(primDiff >= 0)
  timestepLength <- new.magpie(years = getYears(landHarmonized)[-1],
                               fill = diff(getYears(landHarmonized, as.integer = TRUE)))
  stopifnot(timestepLength > 0)
  primHarv <- primDiff / timestepLength
  stopifnot(primHarv >= 0)
  # - 10^-10 to ensure timestepLength * primHarv <= primDiff
  primHarv <- pmax(primHarv - 10^-10, 0)
  stopifnot(timestepLength * primHarv <= primDiff)

  # calculate raw (inconsistent with land) harmonized wood harvest area
  xInput <- calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE)
  xInput <- xInput[, , "wood_harvest_area"]

  xTarget <- calcOutput("NonlandTargetExtrapolated", input = input, target = target,
                        harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
  xTarget <- xTarget[, , "wood_harvest_area"]
  harmonizer <- toolGetHarmonizer(harmonization, constantTotal = FALSE)
  rawHarvestHarmonized <- harmonizer(xInput, xTarget, harmonizationPeriod = harmonizationPeriod)

  stopifnot(rawHarvestHarmonized >= 0)
  if (harmonization != "absoluteChanges") {
    afterHarmonization <- getYears(xInput, as.integer = TRUE)
    afterHarmonization <- afterHarmonization[afterHarmonization >= harmonizationPeriod[2]]
    stopifnot(rawHarvestHarmonized[, afterHarmonization, ] == xInput[, afterHarmonization, ])
  }
  rawHarvestHarmonized <- toolAggregateWoodHarvest(rawHarvestHarmonized)

  histYears <- getYears(xTarget, as.integer = TRUE)
  histYears <- histYears[histYears <= harmonizationPeriod[1]]
  futureYears <- setdiff(getYears(xInput, as.integer = TRUE), histYears)

  if (harmonization == "absoluteChanges") {
    hp <- harmonizationPeriod[1]
    primHarv <- primHarv[, futureYears, ]
    inputRegion <- dimSums(toolAggregateWoodHarvest(xInput), 1.2)
    targetRegion <- dimSums(toolAggregateWoodHarvest(xTarget), 1.2)
    maxHarv <- toolMaxHarvestPerYear(landHarmonized, split = FALSE)[, futureYears, ]
    # the harvest covers, forest and other land, each with its prim source and
    # its secondary/plantation pool, expectation and rescaling work per cover
    covers <- list(forest = list(prim = "primf", secd = c("secdf", "pltns")),
                   other = list(prim = "primn", secd = "secdn"))
    coverExpected <- list()
    coverGlobal <- list()
    coverItems <- list()
    coverTotal <- list()
    coverCap <- list()

    # distribute the cover total of every cell back onto its items, keeping the
    # proportional split of the rescaled amounts as far as the item caps allow,
    # cells that held nothing of the cover fill proportionally to their caps
    fillItems <- function(xAmount, capAmount, total) {
      target <- pmin(total, dimSums(capAmount, 3))
      amountTotal <- dimSums(xAmount, 3)
      scale <- target / amountTotal
      scale[amountTotal == 0] <- 0
      out <- pmin(xAmount * scale, capAmount)
      spare <- pmax(capAmount - out, 0)
      spareTotal <- dimSums(spare, 3)
      gain <- spare * ((target - dimSums(out, 3)) / spareTotal)
      gain[is.na(gain)] <- 0
      # the floating point noise from the proportional fill must not turn into
      # negative harvest area, it crashes downstream consumers of the output
      return(pmax(out + gain, 0))
    }

    for (cover in names(covers)) {
      prim <- covers[[cover]]$prim
      secd <- covers[[cover]]$secd
      inputCover <- dimSums(inputRegion[, , c(prim, secd)], 3)
      targetCover <- dimSums(targetRegion[, , c(prim, secd)], 3)
      expected <- setYears(targetCover[, hp, ], NULL) +
        (inputCover[, futureYears, ] - setYears(inputCover[, hp, ], NULL))
      # global first allocation: floor the regional expectations at zero and
      # scale the rest onto the raw global expectation, keeping the global cover
      # total input parallel even where regional expectations had to be floored
      globalExpected <- dimSums(expected, 1)
      positiveTotal <- dimSums(pmax(expected, 0), 1)
      scale <- pmax(globalExpected, 0) / positiveTotal
      scale[is.na(scale)] <- 0
      coverExpected[[cover]] <- pmax(expected, 0) * scale
      coverGlobal[[cover]] <- globalExpected
    }

    # the negative global cover expectation floored to zero is charged to the
    # other cover, so the global total over all covers stays input parallel
    globalTotal <- coverGlobal$forest + coverGlobal$other
    positiveGlobalTotal <- pmax(coverGlobal$forest, 0) + pmax(coverGlobal$other, 0)
    balanceFactor <- pmax(globalTotal, 0) / positiveGlobalTotal
    balanceFactor[is.na(balanceFactor)] <- 1
    coverGlobalBalanced <- pmax(globalTotal, 0)
    for (cover in names(covers)) {
      coverExpected[[cover]] <- coverExpected[[cover]] * balanceFactor
    }

    for (cover in names(covers)) {
      secd <- covers[[cover]]$secd
      requiredSecd <- coverExpected[[cover]] - dimSums(dimSums(primHarv[, , covers[[cover]]$prim], 1.2), 3)
      requiredSecd[requiredSecd < 0] <- 0

      # the zero clamping of the raw absolute changes and an over-achieving prim
      # harvest broke the parallelism, rescale the cover's secd/pltns totals onto
      # its requirement, keeping the per-cell item split of the raw data
      msgHold <- paste0("Regions with expected ", cover, " harvest area have harvestable land left to hold it")
      secdHarv <- rawHarvestHarmonized[, futureYears, secd]
      secdTotal <- toolRescaleRegionTotals(dimSums(secdHarv, 3), requiredSecd, msgHold)
      share <- secdHarv / dimSums(secdHarv, 3)
      share[is.na(share)] <- 0
      secdHarv <- share * secdTotal

      # move cover harvest above the physical cap into the headroom of the cover
      # within the region, or of the global pool for the years where no region
      # can absorb its own excess
      msgAbsorb <- paste0("Regions absorb their own ", cover, " wood harvest area above the available land")
      coverTotal[[cover]] <- toolRedistributeOverCap(dimSums(secdHarv, 3),
                                                     dimSums(maxHarv[, , secd], 3), msgAbsorb)
      coverCap[[cover]] <- dimSums(maxHarv[, , secd], 3)
      coverItems[[cover]] <- secdHarv
    }

    # harvest of one cover that neither the cover's regional nor its global pool
    # can host is shifted into the headroom of the other cover, region first and
    # then globally, only what remains unhostable anywhere is reported as excess
    for (cover in names(covers)) {
      other <- setdiff(names(covers), cover)
      for (margin in c(1.2, 1)) {
        excess <- pmax(coverTotal[[cover]] - coverCap[[cover]], 0)
        headroom <- pmax(coverCap[[other]] - coverTotal[[other]], 0)
        hostable <- dimSums(headroom, margin) / dimSums(excess, margin)
        hostable[is.na(hostable)] <- 0
        hostable[hostable > 1] <- 1
        filled <- dimSums(excess, margin) / dimSums(headroom, margin)
        filled[is.na(filled)] <- 0
        filled[filled > 1] <- 1
        coverTotal[[cover]] <- coverTotal[[cover]] - excess * hostable
        coverTotal[[other]] <- coverTotal[[other]] + headroom * filled
      }
    }

    # split the final cover totals back onto the items within their caps
    secdHarv <- mbind(fillItems(coverItems$forest, maxHarv[, , covers$forest$secd], coverTotal$forest),
                      fillItems(coverItems$other, maxHarv[, , covers$other$secd], coverTotal$other))
    excessSecdHarvest <- pmax(coverTotal$forest - coverCap$forest, 0) +
      pmax(coverTotal$other - coverCap$other, 0)
  } else {
    # shift harvest from prim to secd to match prim land reduction
    primHarvDiff <- rawHarvestHarmonized[, -1, prim] - primHarv
    secdHarv <- rawHarvestHarmonized[, -1, c("secdf", "secdn", "pltns")]
    secd <- c("secdf", "secdn")
    secdHarv[, , secd] <- (secdHarv[, , secd] + magclass::setNames(primHarvDiff, secd))
    secdHarv[secdHarv < 0] <- 0 # harvest from prim already overachieves harvest target
    stopifnot(primHarv[, , "primf"] + secdHarv[, , "secdf"] + 10^-10 >=
                rawHarvestHarmonized[, -1, "primf"] + rawHarvestHarmonized[, -1, "secdf"],
              primHarv[, , "primn"] + secdHarv[, , "secdn"] + 10^-10 >=
                rawHarvestHarmonized[, -1, "primn"] + rawHarvestHarmonized[, -1, "secdn"])

    # check secd excess harvest, try to shift excess secdf to secdn and vice versa
    maxSecdHarv <- toolMaxHarvestPerYear(landHarmonized, split = FALSE)[, , c("secdf", "secdn", "pltns")]
    excessSecdHarvest <- secdHarv - maxSecdHarv
    excessSecdHarvest[excessSecdHarvest < 0] <- 0
    excessSecdHarvest <- dimSums(excessSecdHarvest, 3)
    secdHarv <- pmin(secdHarv, maxSecdHarv)

    potentialHarvestLeft <- maxSecdHarv - secdHarv
    stopifnot(potentialHarvestLeft >= 0)
    smallest <- pmin(collapseDim(potentialHarvestLeft[, , "secdf"]),
                     collapseDim(potentialHarvestLeft[, , "secdn"]),
                     collapseDim(potentialHarvestLeft[, , "pltns"]))
    stopifnot(smallest[excessSecdHarvest > 0] == 0)

    secdHarv <- pmin(secdHarv + excessSecdHarvest, maxSecdHarv)
    excessSecdHarvest <- excessSecdHarvest - dimSums(potentialHarvestLeft, 3)
    excessSecdHarvest[excessSecdHarvest < 0] <- 0
  }

  # report excess harvest
  excessGlobal <- dimSums(excessSecdHarvest[, futureYears, ], 1)
  excessGlobal <- as.data.frame(excessGlobal, rev = 3)
  excessGlobal <- excessGlobal[excessGlobal$.value > 10^-10, ]
  msg <- paste0(" (global excess harvest area in Mha yr-1: ",
                paste0(excessGlobal$year, ": ", signif(excessGlobal$.value, 3), collapse = ", "),
                "; other years are ok)")
  toolExpectTrue(nrow(excessGlobal) == 0,
                 paste0("Harvest area fits into available land",
                        if (nrow(excessGlobal) != 0) msg))

  # assemble output
  out <- mbind(primHarv, secdHarv)
  out <- out[, futureYears, ]
  map <- toolWoodHarvestMapping()
  map$harvest <- sub("wood_harvest_area\\.", "", map$harvest)
  out <- toolAggregate(out, map,
                       weight = xTarget[, harmonizationPeriod[1], ],
                       from = "land", to = "harvest", dim = 3, zeroWeight = "fix")
  out <- add_dimension(out, dim = 3.1, add = "category", nm = "wood_harvest_area")

  out <- mbind(xTarget[, histYears, ], out)
  if (harmonization != "absoluteChanges") {
    stopifnot(dimSums(out[, -1, ], 3) + excessSecdHarvest + 10^-10 >=
                dimSums(rawHarvestHarmonized[, -1, ], 3))
  }

  # consistency checks
  toolExpectTrue(min(out) >= 0, "All values are >= 0")
  toolExpectLessDiff(out[, histYears, ], xTarget[, histYears, ], 0,
                     "history is not changed by harmonization")
  if (harmonization == "absoluteChanges") {
    outRegionTotal <- dimSums(toolAggregateWoodHarvest(out[, futureYears, ]), 1.2)
    outCoverTotals <- list()
    for (cover in names(covers)) {
      items <- c(covers[[cover]]$prim, covers[[cover]]$secd)
      outCoverTotals[[cover]] <- dimSums(outRegionTotal[, , items], 3)
      message <- paste0("Regional ", cover, " wood harvest area is parallel to input ",
                        "after the harmonization year")
      toolExpectLessDiff(outCoverTotals[[cover]], coverExpected[[cover]], 10^-5, message)
    }
    toolExpectLessDiff(dimSums(outCoverTotals$forest + outCoverTotals$other, 1),
                       coverGlobalBalanced, 10^-5,
                       "Global wood harvest area is parallel to input after the harmonization year")
  }
  toolCheckWoodHarvestArea(out, landHarmonized, harmonizationPeriod[1])

  return(list(x = out,
              isocountries = FALSE,
              unit = "Mha yr-1",
              min = 0,
              description = "Harmonized wood harvest area data"))
}
