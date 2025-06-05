#' calcNonlandHarmonized
#'
#' Harmonize nonland input data to target data using the specified method, checking
#' data for consistency before returning.
#'
#' Wood harvest biomass (bioh) is adapted to the harmonized wood harvest area
#' by calculating kg C per mega hectare for input and target data and
#' harmonizing it. This is then multiplied by the harmonized wood harvest area
#' and scaled so the total harmonized bioh is reached.
#'
#' @param input name of the input dataset, currently only "magpie"
#' @param target name of the target dataset, currently only "luh2"
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param harmonization harmonization method, see \code{\link{toolGetHarmonizer}} for available methods
#' @return harmonized nonland data
#' @author Pascal Sauer
calcNonlandHarmonized <- function(input, target, harmonizationPeriod, harmonization) {
  xInput <- calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE)
  geometry <- attr(xInput, "geometry")
  crs <- attr(xInput, "crs")

  biohMap <- toolBiohMapping()
  kgCPerMhaInput <- xInput[, , biohMap$bioh] / magclass::setNames(xInput[, , woodHarvestAreaCategories()],
                                                                  sub("wood_harvest_area$", "bioh",
                                                                      woodHarvestAreaCategories()))
  kgCPerMhaInput[is.nan(kgCPerMhaInput)] <- 0
  # should we prevent bioh without harvest area before this point, so the following line is no longer necessary?
  kgCPerMhaInput[is.infinite(kgCPerMhaInput)] <- max(kgCPerMhaInput[is.finite(kgCPerMhaInput)])
  stopifnot(is.finite(kgCPerMhaInput), kgCPerMhaInput >= 0)
  getItems(kgCPerMhaInput, 3) <- sub("bioh$", "kgC_per_Mha", getItems(kgCPerMhaInput, 3))

  xTarget <- calcOutput("NonlandTargetExtrapolated", input = input, target = target,
                        harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)

  kgCPerMhaTarget <- xTarget[, , biohMap$bioh] / magclass::setNames(xTarget[, , woodHarvestAreaCategories()],
                                                                    sub("wood_harvest_area$", "bioh",
                                                                        woodHarvestAreaCategories()))
  kgCPerMhaTarget[is.nan(kgCPerMhaTarget)] <- 0
  # should we prevent bioh without harvest area before this point, so the following line is no longer necessary?
  kgCPerMhaTarget[is.infinite(kgCPerMhaTarget)] <- max(kgCPerMhaTarget[is.finite(kgCPerMhaTarget)])
  stopifnot(is.finite(kgCPerMhaTarget), kgCPerMhaTarget >= 0)
  getItems(kgCPerMhaTarget, 3) <- sub("bioh$", "kgC_per_Mha", getItems(kgCPerMhaTarget, 3))

  harmonizer <- toolGetHarmonizer(harmonization)
  out <- harmonizer(mbind(xInput[, , woodHarvestAreaCategories(), invert = TRUE],
                          kgCPerMhaInput),
                    mbind(xTarget[, , woodHarvestAreaCategories(), invert = TRUE],
                          kgCPerMhaTarget),
                    harmonizationPeriod = harmonizationPeriod)

  harvestArea <- calcOutput("WoodHarvestAreaHarmonized", input = input, target = target,
                            harmonizationPeriod = harmonizationPeriod, harmonization = harmonization, aggregate = FALSE)

  # adapt bioh to harmonized harvest area
  kgCPerMhaHarmonized <- out[, , getItems(kgCPerMhaTarget, 3)]
  # should we prevent harvest area without bioh before this point, so the following line is no longer necessary?
  kgCPerMhaHarmonized[kgCPerMhaHarmonized == 0] <- min(kgCPerMhaHarmonized[kgCPerMhaHarmonized > 0])
  stopifnot(is.finite(kgCPerMhaHarmonized), kgCPerMhaHarmonized >= 0)
  biohCalculated <- kgCPerMhaHarmonized * magclass::setNames(harvestArea,
                                                             sub("wood_harvest_area$", "kgC_per_Mha",
                                                                 getItems(harvestArea, 3)))
  getItems(biohCalculated, 3) <- sub("kgC_per_Mha$", "bioh", getItems(biohCalculated, 3))
  stopifnot(0 <= biohCalculated, is.finite(biohCalculated))

  biohNormalization <- dimSums(out[, , biohMap$bioh], 3) / dimSums(biohCalculated, 3)
  biohNormalization[is.nan(biohNormalization)] <- 0
  stopifnot(0 <= biohNormalization, is.finite(biohNormalization))

  biohAdapted <- biohNormalization * biohCalculated

  toolExpectLessDiff(dimSums(out[, , biohMap$bioh], 3),
                     dimSums(biohAdapted, 3),
                     10^-4, "Adapting bioh to harmonized wood harvest area does not change total bioh")

  out[, , biohMap$bioh] <- biohAdapted
  out <- mbind(out, harvestArea)
  out <- out[, , getItems(kgCPerMhaTarget, 3), invert = TRUE]

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
