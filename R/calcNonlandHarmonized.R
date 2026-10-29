#' calcNonlandHarmonized
#'
#' Harmonize nonland input data to target data using the specified method, checking
#' data for consistency before returning.
#'
#' Wood harvest biomass (bioh) is adapted to the harmonized wood harvest area
#' by calculating kg C per mega hectare for input and target data and
#' harmonizing it. This is then multiplied by the harmonized wood harvest area,
#' and the recalculated bioh drives the per-cluster split of the output. For
#' absoluteChanges the regional forest and other land cover totals of bioh stay
#' parallel to the input (constant offset) after the harmonization year: each
#' cover's regional expectation is floored at zero and scaled onto the global
#' expectation, so a region expecting a loss receives nothing while the rest
#' keep the global cover total parallel, and is then allocated onto the
#' recalculated bioh directly (toolAllocateRegionTotals), proportional to
#' harmonized harvest area times harvest rate, so clusters without harvest area
#' host nothing and a region expecting bioh without any harvest area left to
#' hold it hands its share to the global pool, where the receiving regions
#' exceed their own expectation, reported only by the warn-only checks. If the
#' global expectation of a cover itself turns negative, it is floored to zero
#' and the clipped loss is charged to the other covers by uniformly scaling
#' their targets, keeping the global total over all covers parallel. For the
#' other harmonizations bioh of clusters whose harmonized harvest area is zero
#' is redistributed to clusters with harvest area in the same region, falling
#' back to the global pool (toolRedistributeBioh).
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

  if (harmonization == "absoluteChanges") {
    # bioh is not harmonized per cluster, its regional cover expectations are
    # allocated onto the recalculated bioh below instead
    harmonizationInput <- harmonizationInput[, , "bioh", invert = TRUE]
    harmonizationTarget <- harmonizationTarget[, , "bioh", invert = TRUE]
  }

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
    covers <- list(forest = c("bioh.primf", "bioh.secmf", "bioh.secyf", "bioh.pltns"),
                   other = c("bioh.primn", "bioh.secnf"))
    expectedCoverBioh <- list()
    coverGlobal <- list()
    coverItems <- list()
    for (cover in names(covers)) {
      items <- intersect(covers[[cover]], getItems(biohCalculated, 3))
      if (length(items) == 0) {
        next
      }
      inputCoverBioh <- dimSums(xInput[, c(hp[1], yearsAfterHy), items], c(1.2, 3))
      expected <- setYears(dimSums(xTarget[, hp[1], items], c(1.2, 3)), NULL) +
        (inputCoverBioh[, yearsAfterHy, ] - setYears(inputCoverBioh[, hp[1], ], NULL))
      # global first allocation: floor the regional cover expectations at zero
      # and scale the rest onto the raw global expectation, keeping the global
      # cover total parallel even where regional expectations had to be floored
      globalExpected <- dimSums(expected, 1)
      positiveTotal <- dimSums(pmax(expected, 0), 1)
      scale <- pmax(globalExpected, 0) / positiveTotal
      scale[is.na(scale)] <- 0
      expectedCoverBioh[[cover]] <- pmax(expected, 0) * scale
      coverGlobal[[cover]] <- globalExpected
      coverItems[[cover]] <- items
    }

    # the negative global cover expectation floored to zero is charged to the
    # other covers, so the global total over all covers stays input parallel
    globalTotal <- Reduce(`+`, coverGlobal)
    positiveGlobalTotal <- Reduce(`+`, lapply(coverGlobal, function(globalExpected) pmax(globalExpected, 0)))
    balanceFactor <- pmax(globalTotal, 0) / positiveGlobalTotal
    balanceFactor[is.na(balanceFactor)] <- 1
    coverGlobalBalanced <- pmax(globalTotal, 0)
    for (cover in names(expectedCoverBioh)) {
      expectedCoverBioh[[cover]] <- expectedCoverBioh[[cover]] * balanceFactor
    }

    # allocate every cover's regional expectation directly onto the
    # recalculated bioh, proportional to its weight (split onto the items by
    # their recalculated values, which together is proportional to them), so
    # clusters without harvest area host nothing and the regional cover totals
    # hold by construction
    allocatedCoverBioh <- list()
    for (cover in names(expectedCoverBioh)) {
      items <- coverItems[[cover]]
      coverBiohWeight <- dimSums(biohCalculated[, yearsAfterHy, items], 3)
      msgCover <- paste0("Regions with expected ", cover, " bioh have harvest area left to hold it")
      itemShare <- biohCalculated[, yearsAfterHy, items] / coverBiohWeight
      itemShare[is.na(itemShare)] <- 0
      allocatedCoverBioh[[cover]] <-
        toolAllocateRegionTotals(expectedCoverBioh[[cover]], coverBiohWeight, msgCover) * itemShare
    }
    biohAdapted <- do.call(mbind, allocatedCoverBioh)
    histYears <- getYears(xTarget, as.integer = TRUE)
    histYears <- histYears[histYears <= hp[1]]
    out <- mbind(out, mbind(xTarget[, histYears, getItems(biohAdapted, 3)], biohAdapted))
  } else {
    # redistribute bioh from clusters without harvest area
    out[, , "bioh"] <- toolRedistributeBioh(out[, , "bioh"], biohCalculated)

    biohNormalization <- dimSums(out[, , "bioh"], 3) / dimSums(biohCalculated, 3)
    biohNormalization[is.nan(biohNormalization)] <- 0
    stopifnot(0 <= biohNormalization, biohNormalization < Inf)

    biohAdapted <- biohNormalization * biohCalculated

    toolExpectLessDiff(dimSums(out[, , "bioh"], 3),
                       dimSums(biohAdapted, 3),
                       10^-4, "Adapting bioh to harmonized wood harvest area does not change total bioh")

    out[, , "bioh"] <- biohAdapted
  }

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
                       coverGlobalBalanced, 10^-5,
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
