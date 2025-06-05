#' calcLandReport
#'
#' Convert the downscaled land use data to the format required by the given project.
#'
#' @param outputFormat format in which the outputs should be prepared. Options: ESM, ScenarioMIP
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @param yearsSubset vector of years to keep in the output dataset
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @param downscaling name of downscaling method, currently only "magpieClassic"
#' @return land use data
#' @author Pascal Sauer
calcLandReport <- function(outputFormat, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  if (outputFormat == "ESM") {
    native <- calcOutput("LandHighRes", input = "magpie", target = "luh2mod",
                         harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                         downscaling = downscaling, harmonization = harmonization,
                         aggregate = FALSE)
    cellArea <- readSource("LUH2v2h", subtype = "cellArea", convert = FALSE)
    cellArea <- collapseDim(as.magpie(cellArea), 3)

    cropData <- toolCropData(native, cellArea)
    getNames(cropData) <- sub("^cpbf1_", "crpbf_", getNames(cropData))
    getNames(cropData) <- sub("^cpbf2_", "crpbf2_", getNames(cropData))

    totalSecondaryForest <- dimSums(native[, , c("pltns", "secdf")], dim = 3)
    # calculate manaf (managed forest) = pltns share of total secondary forest
    manaf <- native[, , "pltns"] / totalSecondaryForest
    manaf[totalSecondaryForest == 0] <- 0 # replace NAs introduced by 0 / 0
    getNames(manaf) <- "manaf"

    nonCropData <- native[, , c("primf", "primn", "secdf", "secdn", "urban", "pastr", "range")]
    nonCropData[, , "secdf"] <- totalSecondaryForest
    nonCropData <- nonCropData * 10000 / cellArea[getItems(cellArea, 1) %in% getItems(nonCropData, 1), , ]

    out <- mbind(cropData, nonCropData, manaf)

    if (max(out) > 1.0001) {
      toolStatusMessage("warn", paste0("Some shares are > 1 (max: ", max(out), "), setting those to 1"))
    }
    out[out > 1] <- 1

    toolExpectTrue(min(out) >= 0, "All values are >= 0")
    toolExpectTrue(max(out) <= 1, "All shares are <= 1")
    landShares <- c("c3ann", "c3nfx",  "c3per", "c4ann", "c4per", "primf",
                    "primn", "secdf", "secdn", "urban", "pastr", "range")
    toolExpectTrue(max(dimSums(out[, , landShares], 3)) <= 1 + 10^-5,
                   "land shares sum to <= 1 + 10^-5") # if < 1 the rest is likely water
    toolExpectTrue(all(out[, -1, c("primf", "primn")] <= setYears(out[, -nyears(out), c("primf", "primn")],
                                                                  getYears(out[, -1, ]))),
                   "primf and primn are never expanding", falseStatus = "warn")

    return(list(x = out,
                isocountries = FALSE,
                unit = "1",
                min = 0,
                max = 1,
                description = paste("MAgPIE land use data downscaled to LUH2 resolution")))
  } else if (outputFormat == "ScenarioMIP") {
    report <- calcOutput("LandReportScenarioMIP",
                         harmonizationPeriod = harmonizationPeriod,
                         yearsSubset = yearsSubset,
                         harmonization = harmonization, downscaling = downscaling,
                         supplementary = TRUE, aggregate = FALSE)
    return(list(x = report$x,
                isocountries = report$isocountries,
                unit = report$unit,
                min = report$min,
                max = report$max,
                description = report$description))
  } else {
    stop("Can only report for outputFormat = ESM/ScenarioMIP")
  }
}
