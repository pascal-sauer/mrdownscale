#' readMagpieFulldataGdx
#'
#' Read function for data coming from the MAgPIE model.
#'
#' @param subtype type of data to be read in. Available options are
#' land, crop, woodHarvestWeight, woodHarvestArea, fertilizer, clustermap
#' @author Pascal Sauer, Jan Philipp Dietrich
readMagpieFulldataGdx <- function(subtype) {
  "!# @monitor magpie4:::land"
  "!# @monitor magpie4:::croparea"
  "!# @monitor magpie4:::woodProduction"
  "!# @monitor magpie4:::woodHarvestArea"
  "!# @monitor magpie4:::NitrogenBudget"

  gdx <- "fulldata.gdx"
  stopifnot(file.exists(gdx),
            length(Sys.glob("clustermap_*.rds")) == 1)

  clustermap <- readRDS(Sys.glob("clustermap_*.rds"))

  if (subtype == "land") {
    x <- magpie4::land(gdx, level = "cell", subcategories = "crop")
    getSets(x) <- c("region", "id", "year", "data")
    unit <- "Mha"
    description <- "land use information"
  } else if (subtype == "crop") {
    x <- magpie4::croparea(gdx, level = "cell", product_aggr = FALSE, water_aggr = FALSE)
    attr(x, "gdxdata") <- NULL
    getSets(x) <- c("region", "id", "year", "data", "water")
    unit <- "Mha"
    description <- "crop land information separated by irrigated/rainfed"
  } else if (subtype == "woodHarvestWeight") {
    x <- magpie4::woodProduction(gdx)
    stopifnot(identical(getComment(x), " unit: Pg DM yr-1"))
    unit <- "Pg DM yr-1"
    description <- "roundwood and fuelwood harvest weight separated by source"
  } else if (subtype == "woodHarvestArea") {
    x <- magpie4::woodHarvestArea(gdx)
    stopifnot(identical(getComment(x), " unit: Mha yr-1"))
    unit <- "Mha yr-1"
    description <- "wood harvest area separated by source and age classes"
  } else if (subtype == "fertilizer") {
    suppressSpecificWarnings({
      suppressMessages({
        x <- magpie4::NitrogenBudget(gdx, level = "cell", cropTypes = TRUE)
      })
    }, "due to non-iteration of fertilizer distribution, residual fertilizer deficit is moved to balanceflow.")
    x <- collapseDim(x[, , "fertilizer"])
    getSets(x) <- c("region", "id", "year", "cropType")
    unit <- "Tg yr-1"
    description <- "fertilization rate per croptype"
  } else if (subtype == "clustermap") {
    return(list(x = clustermap, class = "data.frame"))
  } else {
    stop("Unknown subtype '", subtype, "' in readMagpieFulldataGdx")
  }

  attr(x, "gdxMetadata") <- NULL

  return(list(x = x, min = 0, unit = unit, description = description))
}
