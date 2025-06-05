#' calcNonlandHighRes
#'
#' Calculate a high resolution dataset from the low resolution input dataset
#' and high resolution data.
#'
#' Wood harvest area is disaggregated using the maximum possible harvest per
#' year which is based on the downscaled land data. Bioh
#' is disaggregated using the just disaggregated wood harvest area as weight.
#' Fertilizer is disaggregated using the downscaled land data. Harvest weight
#' type is disaggregated using the nonland target data in the first year of
#' the harmonization period as weight.
#'
#' @param input name of an input dataset, currently only "magpie"
#' @param target name of a target dataset, currently only "luh2mod"
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method, currently only "magpieClassic"
#' @return downscaled nonland data
#'
#' @examples
#' \dontrun{
#'   calcOutput("NonlandHighRes", input = "magpie", target = "luh2mod",
#'              harmonizationPeriod = c(2015, 2050), yearsSubset = 2015:2100,
#'              harmonization = "fade", downscaling = "magpieClassic")
#' }
#' @author Pascal Sauer
calcNonlandHighRes <- function(input, target, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  x <- calcOutput("NonlandHarmonized", input = input, target = target,
                  harmonizationPeriod = harmonizationPeriod, harmonization = harmonization, aggregate = FALSE)
  x <- x[, getYears(x, as.integer = TRUE) %in% yearsSubset, ]

  futureYears <- getYears(x, as.integer = TRUE)
  futureYears <- futureYears[futureYears > harmonizationPeriod[1]]

  resmap <- calcOutput("ResolutionMapping", input = input, target = target, aggregate = FALSE)

  landHighRes <- calcOutput("LandHighRes", input = input, target = target,
                            harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                            harmonization = harmonization, downscaling = downscaling, aggregate = FALSE)
  land <- landHighRes[, , c("urban", "pastr", "range"), invert = TRUE]
  map <- as.data.frame(rbind(c("primf", "primf"),
                             c("pltns", "pltns"),
                             c("secdf", "secdf"),
                             c("primn", "primn"),
                             c("secdn", "secdn"),
                             c("c3ann_irrigated", "c3ann"),
                             c("c3ann_rainfed", "c3ann"),
                             c("c3ann_irrigated_biofuel_1st_gen", "c3ann"),
                             c("c3ann_rainfed_biofuel_1st_gen", "c3ann"),
                             c("c3nfx_irrigated", "c3nfx"),
                             c("c3nfx_rainfed", "c3nfx"),
                             c("c3nfx_irrigated_biofuel_1st_gen", "c3nfx"),
                             c("c3nfx_rainfed_biofuel_1st_gen", "c3nfx"),
                             c("c3per_irrigated", "c3per"),
                             c("c3per_rainfed", "c3per"),
                             c("c3per_irrigated_biofuel_1st_gen", "c3per"),
                             c("c3per_rainfed_biofuel_1st_gen", "c3per"),
                             c("c3per_irrigated_biofuel_2nd_gen", "c3per"),
                             c("c3per_rainfed_biofuel_2nd_gen", "c3per"),
                             c("c4ann_irrigated", "c4ann"),
                             c("c4ann_rainfed", "c4ann"),
                             c("c4ann_irrigated_biofuel_1st_gen", "c4ann"),
                             c("c4ann_rainfed_biofuel_1st_gen", "c4ann"),
                             c("c4per_irrigated", "c4per"),
                             c("c4per_rainfed", "c4per"),
                             c("c4per_irrigated_biofuel_1st_gen", "c4per"),
                             c("c4per_rainfed_biofuel_1st_gen", "c4per"),
                             c("c4per_irrigated_biofuel_2nd_gen", "c4per"),
                             c("c4per_rainfed_biofuel_2nd_gen", "c4per")))
  colnames(map) <- c("from", "to")
  land <- toolAggregate(land, map, from = "from", to = "to", dim = 3)

  # wood harvest area, use land in previous year (as that is what's harvested) as weight
  whaCat <- woodHarvestAreaCategories()
  harvestArea <- x[, futureYears, whaCat]

  harvestAreaWeight <- toolMaxHarvestPerYear(land)[, futureYears, ] + 10^-30
  harvestAreaDownscaled <- toolAggregate(harvestArea, resmap, weight = harvestAreaWeight,
                                         from = "lowRes", to = "cell", dim = 1)


  bioh <- x[, futureYears, sub("wood_harvest_area$", "bioh", whaCat)]
  weightBioh <- harvestAreaDownscaled + 10^-30
  getItems(weightBioh, 3) <- sub("wood_harvest_area$", "bioh", getItems(weightBioh, 3))
  stopifnot(getItems(weightBioh, 3) == getItems(bioh, 3))
  biohDownscaled <- toolAggregate(bioh, resmap, weight = weightBioh, from = "lowRes", to = "cell", dim = 1)


  fertilizer <- x[, futureYears, grep("fertilizer$", getItems(x, 3))]
  weightFertilizer <- land[, futureYears, sub("_fertilizer$", "", getItems(fertilizer, 3))] + 10^-30
  getItems(weightFertilizer, 3) <- getItems(fertilizer, 3)
  fertilizerDownscaled <- toolAggregate(fertilizer, resmap, weight = weightFertilizer,
                                        from = "lowRes", to = "cell", dim = 1)

  nonlandTarget <- calcOutput("NonlandTarget", target = target, aggregate = FALSE)
  nonlandTarget <- as.magpie(nonlandTarget[[terra::time(nonlandTarget) <= harmonizationPeriod[1] &
                                              terra::time(nonlandTarget) %in% yearsSubset]])

  harvestType <- x[, futureYears, grep("harvest_weight_type$", getItems(x, 3))]
  weightHarvestType <- nonlandTarget[, harmonizationPeriod[1], getItems(harvestType, 3)]
  weightHarvestType <- collapseDim(weightHarvestType) + 10^-30
  harvestTypeDownscaled <- toolAggregate(harvestType, resmap, weight = weightHarvestType,
                                         from = "lowRes", to = "cell", dim = 1)


  out <- mbind(nonlandTarget,
               mbind(harvestAreaDownscaled, biohDownscaled, fertilizerDownscaled, harvestTypeDownscaled))

  inSum <- dimSums(x, dim = 1)
  outSum <- dimSums(out, dim = 1)
  stopifnot(identical(getYears(inSum), getYears(outSum)),
            setequal(getItems(inSum, 3), getItems(outSum, 3)))

  toolExpectLessDiff(inSum, outSum, 0.1,
                     "No significant global sum difference per category before and after downscaling")
  maxdiff <- max(abs(inSum - outSum) / inSum, na.rm = TRUE)
  toolExpectTrue(maxdiff < 10^-5,
                 paste0("Relative global sum difference per category before and after downscaling < 10^-5 ",
                        "(max relative diff: ", signif(maxdiff, 2), ")"))
  toolExpectTrue(setequal(getItems(out, dim = 3), getItems(x, dim = 3)),
                 "Nonland categories remain unchanged")
  toolExpectTrue(min(out) >= 0, "All values are >= 0")

  toolCheckWoodHarvestArea(out[, getYears(landHighRes), whaCat], landHighRes)

  return(list(x = out,
              min = 0,
              isocountries = FALSE,
              unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg yr-1",
              description = "Downscaled nonland data"))
}
