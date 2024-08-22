#' calcNonlandHarmonized
#'
#' Harmonize nonland input data to target data using the specified method, checking
#' data for consistency before returning.
#'
#' @param input name of the input dataset, currently only "magpie"
#' @param target name of the target dataset, currently only "luh2"
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param method harmonization method, see \code{\link{toolGetHarmonizer}} for available methods
#' @return harmonized nonland data
#' @author Pascal Sauer
calcNonlandHarmonized <- function(input = "magpie", target = "luh2mod",
                                  harmonizationPeriod = c(2015, 2050),
                                  method = "fade") {
  xInput <- calcOutput("NonlandHarmonizedCategories", input = input, aggregate = FALSE)
  geometry <- attr(xInput, "geometry")
  crs <- attr(xInput, "crs")

  xTarget <- calcOutput("NonlandLowRes", input = input, target = target, aggregate = FALSE)



  # landLow <- calcOutput("LandTargetExtrapolated", input = input, target = target,
  #                       harmonizationPeriod = harmonizationPeriod, extrapolate = FALSE, aggregate = FALSE)
  # # TODO check
  # harv <- .aggregateWoodHarvest(xTarget[, , endsWith(getItems(xTarget, 3), "wood_harvest_area")])
  # landLow <- landLow[, , getItems(harv, 3)]
  # dif <- landLow - harv
  # where(dif < 0)$true$individual
  # dif["IND.50", "y1995", "primf"] |> print() # could harvest secdf instead
  # dif["IND.50", "y2015", "secdf"] |> print()
  # dif["IND.54", "y2015", "secdf"] |> print()
  # dif["EUR.43", "y2015", "primf"] |> print() # could harvest secdf instead
  # dif["SSA.176", "y2015", "secdf"] |> print() # ?
  # browser()



  if (method == "fade") {
    # extrapolate
    bioh <- grep("bioh$", getItems(xInput, 3), value = TRUE)
    woodHarvestArea <- grep("wood_harvest_area$", getItems(xInput, 3), value = TRUE)
    harvestWeightType <- grep("harvest_weight_type$", getItems(xInput, 3), value = TRUE)
    fertilizer <- grep("fertilizer$", getItems(xInput, 3), value = TRUE)
    stopifnot(setequal(c(bioh, woodHarvestArea, harvestWeightType, fertilizer), getItems(xInput, 3)))

    inputYears <- getYears(xInput, as.integer = TRUE)
    transitionYears <- inputYears[inputYears > harmonizationPeriod[1] & inputYears < harmonizationPeriod[2]]
    woodHarvest <- toolExtrapolate(xTarget[, , c(bioh, woodHarvestArea, harvestWeightType)], transitionYears,
                                   linearModel = FALSE, fallback = "last")
    fertilizer <- toolExtrapolate(xTarget[, , fertilizer], transitionYears)
    xTargetExtrapolated <- mbind(woodHarvest, fertilizer)
    xTargetExtrapolated[xTargetExtrapolated < 0] <- 0
    xTargetExtrapolated <- mbind(xTarget, xTargetExtrapolated)

    harvEx <- .aggregateWoodHarvest(xTargetExtrapolated[, , endsWith(getItems(xTargetExtrapolated, 3), "wood_harvest_area")])
    landEx <- calcOutput("LandTargetExtrapolated", input = input, target = target,
                         harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
    landEx <- landEx[, , getItems(harvEx, 3)]
    dif <- landEx - harvEx
    a <- where(dif < 0)$true$individual
    a[, 'data'] |> table()
    # TODO check

    # harmonize
    out <- toolHarmonizeFade(xInput, xTargetExtrapolated, harmonizationPeriod = harmonizationPeriod)
    stopifnot(setequal(getItems(out, 3), getItems(xInput, 3)))
    
    
    harvHarm <- .aggregateWoodHarvest(out[, , endsWith(getItems(out, 3), "wood_harvest_area")])
    landHarm <- calcOutput("LandHarmonized", input = input, target = target,
                         harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
    landHarm <- landHarm[, , getItems(harvHarm, 3)]
    dif <- landHarm - harvHarm
    dif[dif < 0] |> summary()
    a <- where(dif < 0)$true$individual
    a[, 'data'] |> table()
  } else {
    harmonizer <- toolGetHarmonizer(method)
    out <- harmonizer(xInput, xTarget, harmonizationPeriod = harmonizationPeriod)
  }

  # account for primf/primn expansion being recategorized to secdf/secdn in calcLandHarmonized
  primfixShares <- attr(calcOutput("LandHarmonized", input = input, target = target,
                                   harmonizationPeriod = harmonizationPeriod,
                                   method = method, aggregate = FALSE), "primfixShares")
  if (any(primfixShares < 1)) {
    toolStatusMessage("note", paste("after harmonization primf/primn expansion was replaced",
                                    "by secdf/secdn, adapting wood harvest accordingly"))

    for (category in c("bioh", "wood_harvest_area")) {
      forest <- c("primf", "secyf", "secmf")
      totalBeforeForest <- dimSums(out[, , paste0(forest, "_", category)], 3)

      primForest <- paste0("primf_", category)
      toSecForest <- out[, , primForest] * (1 - primfixShares[, , "primf"])
      out[, , primForest] <- out[, , primForest] - toSecForest
      stopifnot(out[, , primForest] >= 0)
      secYoung <- paste0("secyf_", category)
      secMature <- paste0("secmf_", category)
      youngShare <- out[, , secYoung] / (out[, , secYoung] + out[, , secMature])
      youngShare[is.na(youngShare)] <- 0.5
      stopifnot(all(0 <= youngShare & youngShare <= 1))
      out[, , secYoung] <- out[, , secYoung] + toSecForest * youngShare
      out[, , secMature] <- out[, , secMature] + toSecForest * (1 - youngShare)

      toolExpectLessDiff(totalBeforeForest,
                         dimSums(out[, , paste0(forest, "_", category)], 3),
                         10^-4, paste0("total wood harvest (", category, ") ",
                                       "from forests is not affected by adaptation"))

      # nonforest
      nonforest <- c("primn", "secnf")
      totalBeforeNonforest <- dimSums(out[, , paste0(nonforest, "_", category)], 3)
      primNonforest <- paste0("primn_", category)
      toSecNonforest <- out[, , primNonforest] * (1 - primfixShares[, , "primn"])
      out[, , primNonforest] <- out[, , primNonforest] - toSecNonforest
      stopifnot(out[, , primNonforest] >= 0)
      secNonforest <- paste0("secnf_", category)
      out[, , secNonforest] <- out[, , secNonforest] + toSecNonforest

      toolExpectLessDiff(totalBeforeNonforest,
                         dimSums(out[, , paste0(nonforest, "_", category)], 3),
                         10^-4, paste0("total wood harvest (", category, ") ",
                                       "from non-forests is not affected by adaptation"))
    }
  }

  categories <- grep("wood_harvest_area$", getItems(out, 3), value = TRUE)
  out[, , categories] <- toolWoodHarvestArea(out[, , categories],
                                             calcOutput("LandHarmonized", input = input, target = target,
                                                        harmonizationPeriod = harmonizationPeriod,
                                                        method = method, aggregate = FALSE),
                                             fix = FALSE)

  attr(out, "geometry") <- geometry
  attr(out, "crs")      <- crs

  # checks
  toolExpectTrue(!is.null(attr(out, "geometry")), "Data contains geometry information")
  toolExpectTrue(!is.null(attr(out, "crs")), "Data contains CRS information")
  toolExpectTrue(identical(unname(getSets(out)), c("region", "id", "year", "data")),
                 "Dimensions are named correctly")
  toolExpectTrue(setequal(getItems(out, dim = 3), getItems(xTarget, dim = 3)),
                 "Nonland categories remain unchanged")
  toolExpectTrue(min(out) >= 0, "All values are >= 0")
  # SpatRaster can hold values up to ~10^40 before replacing with Inf, so check we are well below that
  toolExpectTrue(max(out) < 10^30, "All values are < 10^30")

  return(list(x = out,
              isocountries = FALSE,
              unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg yr-1",
              min = 0,
              description = "Harmonized nonland data"))
}
