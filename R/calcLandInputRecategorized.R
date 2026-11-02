#' calcLandInputRecategorized
#'
#' Computes the land input data in target land categories. Splitting of land
#' categories is performed under use of internal land weights reflecting the
#' prevalence of a certain land category in the given area.
#'
#' Mapping from input to target categories is achieved via a merge of a land input
#' mapping to reference categories and a mapping between land target categories and
#' the same reference categories. Thereby a new source or new target can be supported
#' by supplying a map of that new input and/or target to the reference categories.
#'
#' input = "witch": The "rest" category added in calcLandInput is disaggregated into
#' all missing land variables using the reference dataset included in mrdownscale.
#' Thus the following variables are added, but were not at all part of the
#' input scenario data: range, urban, `c4ann_*`, `c3per_*`, `c3nfx_*`
#' This allows the harmonization and downscaling pipeline to continue with a complete
#' dataset, but these variables should not be used or reported (at the very least
#' they have to be checked extensively).
#'
#' @inheritParams calcLandInput
#' @param target name of the land target source to be used
#' @author Jan Philipp Dietrich, Pascal Sauer
calcLandInputRecategorized <- function(input, target) {
  x <- calcOutput("LandInput", input = input, aggregate = FALSE)

  resolutionMapping <- calcOutput("ResolutionMapping", input = input, target = target, aggregate = FALSE)
  x <- x[unique(resolutionMapping$lowRes), , ]

  resolutionMapping$cluster <- resolutionMapping$lowRes
  "!# @monitor magpie4:::addGeometry"
  x <- magpie4::addGeometry(x, resolutionMapping)

  map <- toolLandCategoriesMapping(input, target)
  # get weights for disaggregation to reference categories
  ref <- calcOutput("LandCategorizationWeight", map = map, geometry = attr(x, "geometry"),
                    crs = attr(x, "crs"), aggregate = FALSE)
  y   <- toolAggregate(x, map, dim = 3, from = "dataInput", to = "merge", weight = ref)
  out <- toolAggregate(y, map, dim = 3, from = "merge",     to = "dataOutput")

  # rename primf/secdf/primn/secdn equivalent dataOutput categories to primf/secdf/primn/secdn
  # save original names in woodlandMap to rename back in the end
  woodlandMap <- list()
  for (name in c("primn", "secdn", "primf", "secdf")) {
    woodlandMap[name] <- map$dataOutput[map$reference == name]
    if (sum(map$dataOutput == woodlandMap[name]) == 1) {
      stopifnot(name == woodlandMap[name] || !name %in% getItems(out, 3))
      getItems(out, 3) <- sub(woodlandMap[name], name, getItems(out, 3), fixed = TRUE)
    } else {
      woodlandMap[name] <- NA_character_
    }
  }

  # in some cases other land is mapped to primn and secdn
  # category remapping does not take into account that primn cannot expand, so redistribute
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

  attr(out, "crs") <- attr(x, "crs")
  attr(out, "geometry") <- attr(x, "geometry")

  toolExpectTrue(identical(unname(getSets(out)), c("region", "id", "year", "data")),
                 "Dimensions are named correctly")
  expectedCategories <- setdiff(map$dataOutput, woodlandMap)
  expectedCategories <- c(expectedCategories, names(woodlandMap)[!is.na(woodlandMap)])
  toolExpectTrue(setequal(getItems(out, dim = 3), expectedCategories), "Land categories match target definition")
  toolExpectTrue(all(out >= 0), "All values are >= 0")
  outSum <- dimSums(out, dim = 3)
  toolExpectLessDiff(outSum, outSum[, 1, ], 10^-6, "Total areas stay constant over time")
  toolExpectLessDiff(outSum, dimSums(x, dim = 3), 10^-6,
                     "Total areas are not affected by recategorization")
  toolPrimExpansionCheck(out)

  return(list(x = out,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Input data with land categories remapped to categories of target dataset",
              woodlandMap = woodlandMap,
              clean_magpie = FALSE)) # preserve region ids
}
