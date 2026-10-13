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
#' in kg ha-1 yr-1. The rates are rescaled per crop category and year so that
#' the global fertilizer totals from the harmonization year on keep a constant
#' offset to the input (parallel year-to-year changes), redistributing amounts
#' lost where zero clamping or missing cropland removed fertilizer.
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

  harmonizer <- toolGetHarmonizer(harmonization, constantTotal = FALSE)
  out <- harmonizer(harmonizationInput, harmonizationTarget, harmonizationPeriod = hp)

  landHarmonizedMha <- calcOutput("LandHarmonized", input = input, target = target,
                                  harmonizationPeriod = hp,
                                  harmonization = harmonization, aggregate = FALSE)
  out[, , "fertilizer"] <- toolFertilizerKgPerHa(out[, , "fertilizer"], landHarmonizedMha)

  yearsAfterHy <- getYears(out, as.integer = TRUE)[getYears(out, as.integer = TRUE) > hp[1]]

  if (harmonization == "absoluteChanges") {
    inputGlobal <- dimSums(xInput[, c(hp[1], yearsAfterHy), "fertilizer"], 1)
    expectedGlobalTg <- setYears(dimSums(fertilizerTgTarget[, hp[1], ], 1), NULL) +
      (inputGlobal[, yearsAfterHy, ] - setYears(inputGlobal[, hp[1], ], NULL))
    expectedGlobalTg[expectedGlobalTg < 0] <- 0

    # where the zero clamping of the harmonizer removed all cropland of a cell,
    # the fertilizer rate is undefined, so set it to zero like toolFertilizerKgPerHa
    # does for cells without fertilizer
    fertilizer <- out[, , "fertilizer"]
    fertilizer[fertilizer == Inf] <- 0

    # zero clamping and the Inf -> 0 replacement break the parallelism of the raw
    # absolute changes, rescale the rates per crop category and year so that the
    # global totals equal target total + input delta again, which redistributes
    # fertilizer lost in cells without cropland to the cells that still have some
    actualGlobalTg <- dimSums(toolFertilizerTg(fertilizer[, yearsAfterHy, ], landHarmonizedMha), 1)
    factor <- expectedGlobalTg / actualGlobalTg
    factor[actualGlobalTg == 0] <- 1
    fertilizer[, yearsAfterHy, ] <- fertilizer[, yearsAfterHy, ] * factor
    out[, , "fertilizer"] <- fertilizer
  }

  harvestArea <- calcOutput("WoodHarvestAreaHarmonized", input = input, target = target,
                            harmonizationPeriod = hp, harmonization = harmonization, aggregate = FALSE)

  # adapt bioh to harmonized harvest area
  kgCPerMhaHarmonized <- out[, , getItems(kgCPerMhaTarget, 3)]
  if (harmonization == "absoluteChanges") {
    # the zero clamping of the harmonizer can drive rates to zero, replace them
    # with the smallest positive rate like the floors for input and target data
    kgCPerMhaHarmonized[kgCPerMhaHarmonized == 0] <- min(kgCPerMhaHarmonized[kgCPerMhaHarmonized > 0])
  }
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

  # for years after harmonization make sure the global fertilizer totals match
  # the input, except for absoluteChanges, which deviates from the input by
  # design (target in the harmonization year plus input deltas): there the
  # global totals per crop category must instead stay parallel to the input,
  # the deviation from the input remains as a constant offset
  if (harmonization == "absoluteChanges") {
    fertilizerOutput <- dimSums(toolFertilizerTg(out[, yearsAfterHy, "fertilizer"], landHarmonizedMha), 1)
    toolExpectLessDiff(fertilizerOutput, expectedGlobalTg, 10^-5,
                       "Global fertilizer totals are parallel to input after the harmonization year")
  } else {
    years <- getYears(out, TRUE)[getYears(out, TRUE) >= hp[2]]
    fertilizerInput <- nonlandInput[, years, "fertilizer"]
    fertilizerOutput <- toolFertilizerTg(out[, years, "fertilizer"], landHarmonizedMha[, years, ])
    toolExpectLessDiff(fertilizerInput, fertilizerOutput, 10^-5,
                       "Fertilizer after harmonization period matches input data")
  }
  toolCheckFertilizer(out[, , "fertilizer"])

  return(list(x = out,
              isocountries = FALSE,
              unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg ha-1 yr-1",
              min = 0,
              description = "Harmonized nonland data"))
}
