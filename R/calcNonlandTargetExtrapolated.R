#' calcNonlandTargetExtrapolated
#'
#' Aggregated low resolution target data is extrapolated to the given years
#' using toolExtrapolate. To extrapolate wood harvest weight (bioh) multiply
#' wood harvest area already extrapolated by calcLandTargetExtrapolated with
#' the historical wood harvest rate in kg C per Mha. Fertilizer is extrapolated
#' and returned in kg ha-1 yr-1.
#'
#' @inheritParams calcNonlandInput
#' @param target name of a target dataset, see \code{\link{calcLandTarget}}
#' for available target datasets
#' @param harmonizationPeriod Two integer values, will extrapolate to all years
#' present in input data between harmonization start and end year
#' @return extrapolated nonland target data
#'
#' @examples
#' \dontrun{
#'   calcOutput("NonlandTargetExtrapolated", input = "magpie",
#'              target = "luh3pltns", harmonizationPeriod = c(2025, 2050))
#' }
#' @author Pascal Sauer
calcNonlandTargetExtrapolated <- function(input, target, harmonizationPeriod) {
  hp <- harmonizationPeriod
  xInput <- calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE)
  inputYears <- getYears(xInput, as.integer = TRUE)

  xTarget <- calcOutput("NonlandTargetLowRes", input = input, target = target,
                        endOfHistory = hp[1], aggregate = FALSE)

  transitionYears <- inputYears[inputYears > max(getYears(xTarget, TRUE))
                                & inputYears < hp[2]]

  exFertilizer <- toolExtrapolate(xTarget[, , "fertilizer"], transitionYears)
  exFertilizer[exFertilizer < 0] <- 0
  exFertilizer[exFertilizer > max(xTarget[, , "fertilizer"])] <- max(xTarget[, , "fertilizer"])

  # calculate kg C per Mha in historical period
  histBioh <- dimSums(xTarget[, , "bioh"], 2)
  histHarvestArea <- dimSums(xTarget[, , "wood_harvest_area"], 2)
  kgCPerMha <- 1 / histHarvestArea * collapseDim(histBioh) # order is important here for correct dims
  kgCPerMha[is.nan(kgCPerMha)] <- 0
  stopifnot(is.finite(kgCPerMha))

  # get wood harvest area extrapolation, then apply historical kg C per Mha
  xLand <- calcOutput("LandTargetExtrapolated", input = input, target = target,
                      harmonizationPeriod = hp, aggregate = FALSE, supplementary = TRUE)
  stopifnot(xLand$unit == "Mha",
            !is.null(xLand$woodHarvestArea))
  harvestMha <- xLand$woodHarvestArea[, transitionYears, ]

  stopifnot(setequal(getItems(harvestMha, 3), getItems(kgCPerMha, 3)))
  harvestKgC <- harvestMha * kgCPerMha
  getItems(harvestKgC, 3.1)  <- sub("wood_harvest_area", "bioh", getItems(harvestKgC, 3.1))

  harvestType <- toolExtrapolate(xTarget[, , "harvest_weight_type"], transitionYears)
  harvestType[harvestType < 0] <- 0
  harvestType <- harvestType / dimSums(harvestType, 3)
  harvestType[is.nan(harvestType)] <- 0.5
  stopifnot(is.finite(harvestType))
  harvestType <- harvestType * dimSums(harvestKgC, 3)

  out <- mbind(xTarget,
               mbind(exFertilizer, harvestType, harvestMha, harvestKgC))

  toolExpectLessDiff(dimSums(out[, , "bioh"], 3),
                     dimSums(out[, , "harvest_weight_type"], 3),
                     10^5, "Harvest weight types are consistent")
  toolExpectTrue(min(out) >= 0, "All values are >= 0")
  toolExpectLessDiff(out[, getYears(xTarget), ],
                     xTarget,
                     10^-5, "Returning reference data before harmonization period")
  toolCheckFertilizer(out[, , "fertilizer"])

  return(list(x = out,
              isocountries = FALSE,
              unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg ha-1 yr-1",
              min = 0,
              description = "Extrapolated nonland target data for harmonization"))
}
