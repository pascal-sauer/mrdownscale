calcCropareaMZ <- function(outputFormat, input, harmonizationPeriod, yearsSubset, harmonization, downscaling) {
  x <- calcOutput("LandReport", outputFormat = outputFormat, input = input,
                  harmonizationPeriod = harmonizationPeriod, yearsSubset = yearsSubset,
                  harmonization = harmonization, downscaling = downscaling,
                  aggregate = FALSE)

  toolExpectTrue(dim(x)[1] == 67420, "expected 67420 grid cells are present")
  toolExpectTrue(identical(getItems(x, 2),
                           paste0("y", c(seq(1995, 2060, 5), seq(2070, 2100, 10)))),
                 "expected years from 1995 to 2100 are present")

  crops <- c("tece", "maiz", "trce", "rice_pro", "foddr", "soybean", "rapeseed", "groundnut", "sunflower", "oilpalm",
             "puls_pro", "potato", "cassav_sp", "sugr_cane", "sugr_beet", "others", "cottn_pro", "begr", "betr")
  # TODO this check currently fails
  toolExpectTrue(setequal(getItems(x, 3),
                          c(paste0(crops, ".rainfed"), paste0(crops, ".irrigated"))),
                 "expected crops (rainfed and irrigated) are present")

  return(list(x = x,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              max = 1,
              description = paste0("MAgPIE croparea harmonized and downscaled using landuseinit as reference ",
                                   "for further processing in magpie4")))
}
