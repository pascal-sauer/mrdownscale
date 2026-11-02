#' toolCropData
#'
#' For each crop type (c3ann, c3nfx, c3per, c4ann, c4per) calculate irrigation
#' share, 1st and 2nd generation biofuel shares and share of cell area
#' (corresponding to LUH3 variables c3ann, irrig_c3ann, cpbf1_c3ann,
#' cpbf2_c3ann and analogously for the other crop types).
#'
#' @param landHighRes high resolution land use data as magclass object
#' @param cellArea corresponding magclass object containing cell area in Mha
#' @return crop data as magclass object with variables `irrig_*`, `cpbf1_*`,
#' `cpbf2_*`, `*` for all 5 crop types
#' @author Pascal Sauer
toolCropData <- function(landHighRes, cellArea) {
  cropTypes <- c("c3ann", "c3nfx", "c3per", "c4ann", "c4per")
  return(do.call(mbind, lapply(cropTypes, function(cropType) {
    x <- landHighRes[, , grep(cropType, getItems(landHighRes, 3))]
    total <- dimSums(x, dim = 3)

    irrigationShare <- dimSums(x[, , grep("irrigated", getItems(x, 3))], dim = 3) / total
    getNames(irrigationShare) <- paste0("irrig_", cropType)
    irrigationShare[total == 0] <- 0 # replace NAs introduced by 0 / 0

    biofuel1stGenShare <- dimSums(x[, , grep("biofuel_1st_gen", getItems(x, 3))], dim = 3) / total
    getNames(biofuel1stGenShare) <- paste0("cpbf1_", cropType)
    biofuel1stGenShare[total == 0] <- 0 # replace NAs introduced by 0 / 0

    # multiply by 10000 to convert from Mha to km2, divide by cellArea to get shares
    cellAreaShare <- total * 10000 / cellArea[getItems(cellArea, 1) %in% getItems(x, 1), , ]
    getNames(cellAreaShare) <- cropType

    combined <- mbind(irrigationShare, biofuel1stGenShare, cellAreaShare)
    if (any(grepl("biofuel_2nd_gen", getItems(x, 3)))) {
      biofuel2ndGenShare <- dimSums(x[, , grep("biofuel_2nd_gen", getItems(x, 3))], dim = 3) / total
      getNames(biofuel2ndGenShare) <- paste0("cpbf2_", cropType)
      biofuel2ndGenShare[total == 0] <- 0 # replace NAs introduced by 0 / 0
      combined <- mbind(combined, biofuel2ndGenShare)
    }
    return(combined)
  })))
}
