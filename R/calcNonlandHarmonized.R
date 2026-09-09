#' calcNonlandHarmonized
#'
#' Harmonize nonland input data to target data using the specified method, checking
#' data for consistency before returning.
#'
#' Wood harvest biomass (bioh) is adapted to the harmonized wood harvest area
#' by calculating kg C per mega hectare for input and target data and
#' harmonizing it. This is then multiplied by the harmonized wood harvest area
#' and scaled so the total harmonized bioh is reached.
#' Harmonize absolute fertilizer in Tg yr-1, then convert to fertilizer rate
#' in kg ha-1 yr-1.
#'
#' @inheritParams calcNonlandInput
#' @param target name of the target dataset
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param harmonization harmonization method, see \code{\link{toolGetHarmonizer}} for available methods
#' @return harmonized nonland data
#' @author Pascal Sauer
calcNonlandHarmonized <- function(input, target, harmonizationPeriod, harmonization) {
  if (harmonization == "absoluteChanges") {
    stop("Unsupported harmonization method \"", harmonization, "\", ",
         "absoluteChanges is only supported for land harmonization")
  }
  hp <- harmonizationPeriod
  xInput <- calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE)
  geometry <- attr(xInput, "geometry")
  crs <- attr(xInput, "crs")

  kgCPerMhaInput <- xInput[, , "bioh"] / collapseDim(xInput[, , "wood_harvest_area"])
  kgCPerMhaInput[is.nan(kgCPerMhaInput)] <- min(kgCPerMhaInput[!is.nan(kgCPerMhaInput)])
  stopifnot(0 < kgCPerMhaInput, kgCPerMhaInput < Inf)
  getItems(kgCPerMhaInput, 3.1) <- sub("bioh$", "kgC_per_Mha", getItems(kgCPerMhaInput, 3.1))

  xTarget <- calcOutput("NonlandTargetExtrapolated", input = input, target = target,
                        harmonizationPeriod = hp, aggregate = FALSE)

  kgCPerMhaTarget <- xTarget[, , "bioh"] / collapseDim(xTarget[, , "wood_harvest_area"])
  kgCPerMhaTarget[!is.finite(kgCPerMhaTarget)] <- min(kgCPerMhaTarget[is.finite(kgCPerMhaTarget)])
  stopifnot(0 < kgCPerMhaTarget, kgCPerMhaTarget < Inf)
  getItems(kgCPerMhaTarget, 3.1) <- sub("bioh$", "kgC_per_Mha", getItems(kgCPerMhaTarget, 3.1))

  harmonizationInput <- mbind(xInput[, , "wood_harvest_area", invert = TRUE], kgCPerMhaInput)
  fertilizerTgTarget <- toolFertilizerTg(xTarget[, , "fertilizer"],
                                         calcOutput("LandTargetExtrapolated", input = input, target = target,
                                                    harmonizationPeriod = hp, aggregate = FALSE))
  harmonizationTarget <- mbind(xTarget[, , c("wood_harvest_area", "fertilizer"), invert = TRUE],
                               kgCPerMhaTarget, fertilizerTgTarget)

  harmonizer <- toolGetHarmonizer(harmonization)
  out <- harmonizer(harmonizationInput, harmonizationTarget, harmonizationPeriod = hp)

  landHarmonizedMha <- calcOutput("LandHarmonized", input = input, target = target,
                                  harmonizationPeriod = hp,
                                  harmonization = harmonization, aggregate = FALSE)
  out[, , "fertilizer"] <- toolFertilizerKgPerHa(out[, , "fertilizer"], landHarmonizedMha)

  harvestArea <- calcOutput("WoodHarvestAreaHarmonized", input = input, target = target,
                            harmonizationPeriod = hp, harmonization = harmonization, aggregate = FALSE)

  # adapt bioh to harmonized harvest area
  kgCPerMhaHarmonized <- out[, , getItems(kgCPerMhaTarget, 3)]
  stopifnot(0 < kgCPerMhaHarmonized, kgCPerMhaHarmonized < Inf)
  biohCalculated <- kgCPerMhaHarmonized * collapseDim(harvestArea)
  getItems(biohCalculated, 3.1) <- sub("kgC_per_Mha$", "bioh", getItems(biohCalculated, 3.1))
  stopifnot(0 <= biohCalculated, biohCalculated < Inf)

  biohNormalization <- dimSums(out[, , "bioh"], 3) / dimSums(biohCalculated, 3)
  biohNormalization[is.nan(biohNormalization)] <- 0
  stopifnot(0 <= biohNormalization, biohNormalization < Inf)

  biohAdapted <- biohNormalization * biohCalculated

  toolExpectLessDiff(dimSums(out[, , "bioh"], 3),
                     dimSums(biohAdapted, 3),
                     10^-4, "Adapting bioh to harmonized wood harvest area does not change total bioh")

  out[, , "bioh"] <- biohAdapted
  out <- mbind(out, harvestArea)
  out <- out[, , getItems(kgCPerMhaTarget, 3), invert = TRUE]

  attr(out, "geometry") <- geometry
  attr(out, "crs")      <- crs

  nonlandInput <- calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE)

  # checks
  toolExpectTrue(!is.null(attr(out, "geometry")), "Data contains geometry information")
  toolExpectTrue(!is.null(attr(out, "crs")), "Data contains CRS information")
  toolExpectTrue(identical(unname(getSets(out)), c("region", "id", "year", "category", "data")),
                 "Dimensions are named correctly")
  toolExpectTrue(setequal(getItems(out, dim = 3), getItems(xTarget, dim = 3)),
                 "Nonland categories remain unchanged")
  toolExpectTrue(min(out) >= 0, "All values are >= 0")
  # SpatRaster can hold values up to ~10^40 before replacing with Inf, so check we are well below that
  toolExpectTrue(max(out) < 10^30, "All values are < 10^30")

  toolExpectLessDiff(out[, getYears(out, as.integer = TRUE) <= hp[1], ],
                     xTarget[, getYears(xTarget, as.integer = TRUE) <= hp[1], ],
                     10^-4, "Returning reference data before harmonization period")

  # for years after harmonization make sure that total global fertilizer applied matches input
  years <- getYears(out, TRUE)[getYears(out, TRUE) >= hp[2]]
  fertilizerInput <- nonlandInput[, years, "fertilizer"]
  fertilizerOutput <- toolFertilizerTg(out[, years, "fertilizer"], landHarmonizedMha[, years, ])
  toolExpectLessDiff(fertilizerInput, fertilizerOutput, 10^-5,
                     "Fertilizer after harmonization period matches input data")
  toolCheckFertilizer(out[, , "fertilizer"])

  return(list(x = out,
              isocountries = FALSE,
              unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg ha-1 yr-1",
              min = 0,
              description = "Harmonized nonland data"))
}
