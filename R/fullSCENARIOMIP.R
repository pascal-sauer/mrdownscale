fullSCENARIOMIP <- function(rev = numeric_version("0"), ..., scenario = "",
                            harmonizationPeriod = c(2020, 2050),
                            yearsSubset = 2020:2100,
                            harmonization = "fade", downscaling = "magpieClassic",
                            compression = 2, progress = TRUE) {
  stopifnot(...length() == 0)

  revision <- if (identical(rev, numeric_version("0"))) format(Sys.time(), "%Y-%m-%d") else rev

  fileSuffix <- paste0("_input4MIPs_landState_ScenarioMIP_PIK-MAgPIE-",
                       scenario, if (scenario != "") "-",
                       revision, "_gn_", min(yearsSubset), "-", max(yearsSubset), ".nc")

  writeArgs <- list(compression = compression, missval = 1e20, progress = progress,
                    gridDefinition = c(-179.875, 179.875, -89.875, 89.875, 0.25))

  metadataArgs <- list(revision = revision, missingValue = 1e20, resolution = 0.25,
                       compression = compression, harmonizationPeriod = harmonizationPeriod,
                       activityId = "ScenarioMIP",
                       references = paste0("https://github.com/pik-piam/mrdownscale and ",
                                           "https://wcrp-cmip.org/mips/scenariomip/"),
                       targetMIP = "ScenarioMIP",
                       ncTitle = "MAgPIE Land-Use Data Harmonized and Downscaled using LUH3 historic as reference",
                       referenceDataset = paste0("LUH3 historic from https://aims2.llnl.gov/search/input4mips/ ",
                                                 "(institution_id = 'UofMD' and mip_era = 'CMIP7')"),
                       furtherInfoUrl = "NA")

  ncFile <- paste0("multiple-states", fileSuffix)
  calcOutput("StatesNC", outputFormat = "ScenarioMIP", harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  ncFile <- paste0("multiple-management", fileSuffix)
  calcOutput("ManagementNC", outputFormat = "ScenarioMIP",
             harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  ncFile <- paste0("multiple-transitions", fileSuffix)
  calcOutput("TransitionsNC", outputFormat = "ScenarioMIP",
             harmonizationPeriod = harmonizationPeriod,
             yearsSubset = yearsSubset,
             harmonization = harmonization, downscaling = downscaling,
             aggregate = FALSE, file = ncFile, writeArgs = writeArgs)
  do.call(toolAddMetadataNC, c(ncFile = ncFile, metadataArgs))

  toolWriteMadratLog(logPath = "consistencyCheck.log")
}
