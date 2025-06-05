#' calcWoodHarvestAreaHarmonized
#'
#' Harmonize wood harvest area based on harmonized land data.
#' First, wood harvest area is harmonized just like land data.
#' Then, to ensure consistency, primary harvest is converted to seondary
#' harvest (or vice versa) to match the primary land reduction from the
#' harmonized land data. If secondary harvest exceeds the available
#' secondary land, the excess is shifted to the other secondary land
#' (secdf to secdn and vice versa). Remaining excess harvest area is reported.
#'
#' @param input name of the input dataset, currently only "magpie"
#' @param target name of the target dataset, currently only "luh2mod"
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param harmonization harmonization method, see \code{\link{toolGetHarmonizer}} for available methods
#' @return harmonized wood harvest area data
#'
#' @examples
#' \dontrun{
#'   calcOutput("WoodHarvestAreaHarmonized", input = "magpie",
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
  xInput <- xInput[, , woodHarvestAreaCategories()]

  xTarget <- calcOutput("NonlandTargetExtrapolated", input = input, target = target,
                        harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
  xTarget <- xTarget[, , woodHarvestAreaCategories()]
  harmonizer <- toolGetHarmonizer(harmonization)
  rawHarvestHarmonized <- harmonizer(xInput, xTarget, harmonizationPeriod = harmonizationPeriod)

  afterHarmonization <- getYears(xInput, as.integer = TRUE)
  afterHarmonization <- afterHarmonization[afterHarmonization >= harmonizationPeriod[2]]
  stopifnot(rawHarvestHarmonized >= 0,
            rawHarvestHarmonized[, afterHarmonization, ] == xInput[, afterHarmonization, ])
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
  secdHarv <- mpmin(secdHarv, maxSecdHarv)

  potentialHarvestLeft <- maxSecdHarv - secdHarv
  stopifnot(potentialHarvestLeft >= 0)
  smallest <- mpmin(mpmin(potentialHarvestLeft[, , "secdf"],
                          potentialHarvestLeft[, , "secdn"]),
                    potentialHarvestLeft[, , "pltns"])
  stopifnot(smallest[excessSecdHarvest > 0] == 0)

  secdHarv <- mpmin(secdHarv + excessSecdHarvest, maxSecdHarv)
  excessSecdHarvest <- excessSecdHarvest - dimSums(potentialHarvestLeft, 3)
  excessSecdHarvest[excessSecdHarvest < 0] <- 0

  # report excess harvest
  histYears <- getYears(xTarget, as.integer = TRUE)
  histYears <- histYears[histYears <= harmonizationPeriod[1]]
  futureYears <- setdiff(getYears(xInput, as.integer = TRUE), histYears)

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
  out <- toolAggregate(out, toolWoodHarvestMapping(),
                       weight = xTarget[, harmonizationPeriod[1], ] + 10^-30,
                       from = "land", to = "harvest", dim = 3)

  out <- mbind(xTarget[, histYears, ], out)
  stopifnot(dimSums(out[, -1, ], 3) + excessSecdHarvest + 10^-10 >=
              dimSums(rawHarvestHarmonized[, -1, ], 3))

  # consistency checks
  toolExpectTrue(min(out) >= 0, "All values are >= 0")
  toolExpectLessDiff(out[, histYears, ], xTarget[, histYears, ], 0,
                     "history is not changed by harmonization")
  toolCheckWoodHarvestArea(out[, histYears, ], landHarmonized[, histYears, ], "In historical period, ")
  toolCheckWoodHarvestArea(out[, futureYears, ], landHarmonized[, futureYears, ], "After historical period, ")

  return(list(x = out,
              isocountries = FALSE,
              unit = "Mha yr-1",
              min = 0,
              description = "Harmonized wood harvest area data"))
}
