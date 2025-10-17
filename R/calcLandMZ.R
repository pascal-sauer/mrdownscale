calcLandMZ <- function(outputFormat, input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  x <- calcOutput("LandReport", outputFormat = outputFormat, input = input,
                  harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                  harmonization = harmonization, downscaling = downscaling,
                  aggregate = FALSE)

  toolExpectTrue(dim(x)[1] == 67420, "expected 67420 grid cells are present")
  toolExpectTrue(identical(getItems(x, 2),
                           paste0("y", c(1985, seq(1995, 2060, 5), seq(2070, 2100, 10)))),
                 "expected years from 1985 to 2100 are present")

  toolExpectTrue(setequal(getItems(x, 3),
                          c("crop", "past", "forestry", "primforest", "secdforest", "urban", "other")),
                 "expected 7 land types are present")

  return(list(x = x,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              max = 1,
              description = paste0("MAgPIE land data harmonized and downscaled using landuseinit as reference ",
                                   "for further processing in magpie4")))
}
