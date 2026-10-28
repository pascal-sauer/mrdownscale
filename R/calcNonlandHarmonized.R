#' calcNonlandHarmonized
#'
#' Harmonize nonland input data to target data using the specified method, checking
#' data for consistency before returning.
#'
#' Wood harvest biomass (bioh) is adapted to the harmonized wood harvest area
#' by calculating kg C per mega hectare for input and target data and
#' harmonizing it. This is then multiplied by the harmonized wood harvest area
#' and scaled so the total harmonized bioh is reached. For absoluteChanges the
#' regional forest and other land cover totals of bioh stay parallel to the
#' input (constant offset) after the harmonization year: where zero clamping
#' inflated a region's total above its expectation it is rescaled back, the
#' per-cell kg C per mega hectare times harvest area split only deciding which
#' cells can host how much. Because that per-cell split moves bioh between the
#' covers, a final cover-granular rescale restores the parallel total of each
#' cover, falling back to a global factor where a region expects bioh of a
#' cover but has no cells left to hold it. bioh of clusters
#' whose harmonized harvest area is zero is redistributed to clusters with
#' harvest area in the same region, falling back to the global pool
#' (toolRedistributeBioh).
#' Harmonize absolute fertilizer in Tg yr-1, then convert to fertilizer rate
#' in kg ha-1 yr-1. For absoluteChanges the tonnage changes of the input since
#' the harmonization year are applied regionally, projected onto the harmonized
#' crop group areas (toolProjectFertilizerTg). The rates are then rescaled per
#' region, crop category and year so that each region's fertilizer total from
#' the harmonization year on keeps a constant offset to its own input (parallel
#' year-to-year changes), redistributing amounts lost where zero clamping or
#' missing cropland removed fertilizer to the other cells of the same region,
#' falling back to the global factor exactly in those crop group and year
#' tables where a region lost all its cropland and can no longer hold its
#' tonnes, and finally capped at the plausibility threshold
#' (toolCapFertilizerKgPerHa).
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

  yearsAfterHy <- getYears(out, as.integer = TRUE)[getYears(out, as.integer = TRUE) > hp[1]]

  out[, , "fertilizer"] <- toolFertilizerKgPerHa(out[, , "fertilizer"], landHarmonizedMha)

  if (harmonization == "absoluteChanges") {
    # treat the harmonized land as truth for fertilizer tonnages: instead of the
    # per-cell input changes, the regional tonnage changes are projected onto
    # the harmonized crop group areas, so tonnes only grow where the harmonized
    # land still has cropland of that crop group
    fertilizerTgInput <- xInput[, c(hp[1], yearsAfterHy), "fertilizer"]
    projectedTg <- toolProjectFertilizerTg(fertilizerTgTarget, fertilizerTgInput,
                                           landHarmonizedMha, hp[1])
    regionTgInput <- dimSums(fertilizerTgInput, 1.2)
    expectedRegionTg <- setYears(dimSums(fertilizerTgTarget[, hp[1], ], 1.2), NULL) +
      (regionTgInput[, yearsAfterHy, ] - setYears(regionTgInput[, hp[1], ], NULL))
    expectedRegionTg[expectedRegionTg < 0] <- 0
    # the global check below stays exact: it compares against the regional
    # expectations summed up
    expectedGlobalTg <- dimSums(expectedRegionTg, 1)

    # cells whose harmonized land has no cropland of the crop group cannot hold
    # tonnes, they lose them, which the rescaling compensates
    cropMha <- toolAggregateCropland(landHarmonizedMha[getItems(fertilizerTgInput, 1), yearsAfterHy, ],
                                     getItems(fertilizerTgInput, 3.2), keepOthers = FALSE)
    projectedTg[cropMha == 0] <- 0

    # zero clamping and the lost tonnes broke the parallelism, rescale per
    # region, crop category and year to restore each region's expected tonnage
    rescaledTg <- toolRescaleRegionTotals(projectedTg, expectedRegionTg,
                                          "Regions with expected fertilizer tonnage have cropland left to hold it")

    # cap the rates at the plausibility threshold of toolCheckFertilizer, which
    # moves the anchor tonnage left on collapsed-area cells into the headroom
    # of the other cells of their region
    rescaledKgPerHa <- toolFertilizerKgPerHa(rescaledTg, landHarmonizedMha)
    out[, yearsAfterHy, "fertilizer"] <- toolCapFertilizerKgPerHa(rescaledKgPerHa, landHarmonizedMha)
  }

  harvestArea <- calcOutput("WoodHarvestAreaHarmonized", input = input, target = target,
                            harmonizationPeriod = hp, harmonization = harmonization, aggregate = FALSE)

  # adapt bioh to harmonized harvest area
  kgCPerMhaHarmonized <- out[, , getItems(kgCPerMhaTarget, 3)]
  if (harmonization == "absoluteChanges") {
    # floor clamped zero rates to the smallest positive rate, as done for input and target
    kgCPerMhaHarmonized[kgCPerMhaHarmonized == 0] <- min(kgCPerMhaHarmonized[kgCPerMhaHarmonized > 0])
  }
  stopifnot(0 < kgCPerMhaHarmonized, kgCPerMhaHarmonized < Inf)
  biohCalculated <- kgCPerMhaHarmonized * collapseDim(harvestArea)
  getItems(biohCalculated, 3.1) <- sub("kgC_per_Mha$", "bioh", getItems(biohCalculated, 3.1))
  stopifnot(0 <= biohCalculated, biohCalculated < Inf)

  if (harmonization == "absoluteChanges") {
    # zero clamping of the raw absolute changes broke the parallelism, rescale
    # each region's bioh total onto its expectation, keeping the per-cell split
    regionBiohInput <- dimSums(xInput[, c(hp[1], yearsAfterHy), "bioh"], c(1.2, 3))
    expectedRegionBioh <- setYears(dimSums(xTarget[, hp[1], "bioh"], c(1.2, 3)), NULL) +
      (regionBiohInput[, yearsAfterHy, ] - setYears(regionBiohInput[, hp[1], ], NULL))
    expectedRegionBioh[expectedRegionBioh < 0] <- 0
    bioh <- out[, yearsAfterHy, "bioh"]
    biohTotal <- dimSums(bioh, 3)
    rescaledBioh <- toolRescaleRegionTotals(biohTotal, expectedRegionBioh,
                                            "Regions with expected bioh have harvest area left to hold it")
    # the final item split is driven by kgC per Mha times harvest area, so only
    # the totals can be parallel, split the rescaled total back proportionally
    share <- bioh / biohTotal
    share[is.na(share)] <- 0
    out[, yearsAfterHy, "bioh"] <- share * rescaledBioh
  }

  # redistribute bioh from clusters without harvest area
  out[, , "bioh"] <- toolRedistributeBioh(out[, , "bioh"], biohCalculated)

  biohNormalization <- dimSums(out[, , "bioh"], 3) / dimSums(biohCalculated, 3)
  biohNormalization[is.nan(biohNormalization)] <- 0
  stopifnot(0 <= biohNormalization, biohNormalization < Inf)

  biohAdapted <- biohNormalization * biohCalculated

  toolExpectLessDiff(dimSums(out[, , "bioh"], 3),
                     dimSums(biohAdapted, 3),
                     10^-4, "Adapting bioh to harmonized wood harvest area does not change total bioh")

  if (harmonization == "absoluteChanges") {
    # the renormalization re-split the bioh of every cell by harvest area,
    # moving it between the forest and other land covers, restore each region's
    # cover parallelism, keeping the per-cell item split of the renormalized data
    covers <- list(forest = c("bioh.primf", "bioh.secmf", "bioh.secyf", "bioh.pltns"),
                   other = c("bioh.primn", "bioh.secnf"))
    expectedCoverBioh <- list()
    for (cover in names(covers)) {
      items <- intersect(covers[[cover]], getItems(biohAdapted, 3))
      if (length(items) == 0) {
        next
      }
      inputCoverBioh <- dimSums(xInput[, c(hp[1], yearsAfterHy), items], c(1.2, 3))
      expected <- setYears(dimSums(xTarget[, hp[1], items], c(1.2, 3)), NULL) +
        (inputCoverBioh[, yearsAfterHy, ] - setYears(inputCoverBioh[, hp[1], ], NULL))
      expected[expected < 0] <- 0
      expectedCoverBioh[[cover]] <- expected
      coverBioh <- biohAdapted[, yearsAfterHy, items]
      msgCover <- paste0("Regions with expected ", cover, " bioh have harvest area left to hold it")
      coverTotal <- toolRescaleRegionTotals(dimSums(coverBioh, 3), expected, msgCover)
      share <- coverBioh / dimSums(coverBioh, 3)
      share[is.na(share)] <- 0
      biohAdapted[, yearsAfterHy, items] <- share * coverTotal
    }
  }

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

  # absoluteChanges must stay parallel to the input (constant offset) instead of matching it
  if (harmonization == "absoluteChanges") {
    toolExpectLessDiff(dimSums(toolFertilizerTg(out[, yearsAfterHy, "fertilizer"], landHarmonizedMha), 1.2),
                       expectedRegionTg, 10^-5,
                       "Regional fertilizer totals are parallel to input after the harmonization year")
    fertilizerOutput <- dimSums(toolFertilizerTg(out[, yearsAfterHy, "fertilizer"], landHarmonizedMha), 1)
    toolExpectLessDiff(fertilizerOutput, expectedGlobalTg, 10^-5,
                       "Global fertilizer totals are parallel to input after the harmonization year")
    outCoverTotals <- list()
    for (cover in names(covers)) {
      items <- intersect(covers[[cover]], getItems(out, 3))
      if (length(items) == 0) {
        next
      }
      outCoverTotals[[cover]] <- dimSums(out[, yearsAfterHy, items], c(1.2, 3))
      toolExpectLessDiff(outCoverTotals[[cover]], expectedCoverBioh[[cover]], 10^-5,
                         paste0("Regional ", cover, " bioh totals are parallel to input ",
                                "after the harmonization year"))
    }
    toolExpectLessDiff(dimSums(Reduce(`+`, outCoverTotals), 1),
                       dimSums(Reduce(`+`, expectedCoverBioh), 1), 10^-5,
                       "Global bioh totals are parallel to input after the harmonization year")
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
