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
#' harvest area stay parallel to the input (constant offset) and are allocated
#' onto the harvestable land (toolParallelHarvestCoverTotals).
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

  histYears <- getYears(xTarget, as.integer = TRUE)
  histYears <- histYears[histYears <= harmonizationPeriod[1]]
  futureYears <- setdiff(getYears(xInput, as.integer = TRUE), histYears)

  if (harmonization == "absoluteChanges") {
    # the cover targets stay parallel to the input and are allocated onto
    # the harvestable land directly, the checks below use what the tool returns
    primHarv <- primHarv[, futureYears, ]
    coverHarvest <- toolParallelHarvestCoverTotals(xInput, xTarget, primHarv, landHarmonized,
                                                   harmonizationPeriod[1])
    secdHarv <- coverHarvest$secdHarv
    excessSecdHarvest <- coverHarvest$excessSecdHarvest
    coverTarget <- coverHarvest$coverTarget
    coverGlobalBalanced <- coverHarvest$coverGlobalBalanced
    covers <- coverHarvest$covers
  } else {
    harmonizer <- toolGetHarmonizer(harmonization, constantTotal = FALSE)
    rawHarvestHarmonized <- harmonizer(xInput, xTarget, harmonizationPeriod = harmonizationPeriod)
    stopifnot(rawHarvestHarmonized >= 0)
    afterHarmonization <- getYears(xInput, as.integer = TRUE)
    afterHarmonization <- afterHarmonization[afterHarmonization >= harmonizationPeriod[2]]
    stopifnot(rawHarvestHarmonized[, afterHarmonization, ] == xInput[, afterHarmonization, ])
    rawHarvestHarmonized <- toolAggregateWoodHarvest(rawHarvestHarmonized)

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
      msg <- paste0("Regional ", cover, " wood harvest area is parallel to input ",
                    "after the harmonization year")
      toolExpectLessDiff(outCoverTotals[[cover]], coverTarget[[cover]], 10^-5, msg)
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
