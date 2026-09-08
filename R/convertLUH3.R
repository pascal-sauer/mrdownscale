#' convertLUH3
#'
#' Convert LUH3 cell area shares to absolute areas by multiplying with cell area, convert to Mha.
#'
#' @param x SpatRaster with LUH3 cell area shares
#' @param subtype Only "states" is converted
convertLUH3 <- function(x, subtype) {
  if (subtype != "states") {
    stop("for subtype != states pass convert = FALSE")
  }

  cellAreaKm2 <- readSource("LUH3", subtype = "cellArea", subset = NULL, convert = FALSE)
  stopifnot(terra::units(cellAreaKm2) == "km2")
  # convert from km2 to Mha
  cellAreaMha <- cellAreaKm2 / 10000

  stopifnot(max(terra::minmax(x, compute = TRUE)) <= 1.0001,
            all(terra::units(x) == "1"))
  x <- x * cellAreaMha
  unit <- "Mha"
  terra::units(x) <- unit

  # write to disk to save memory, also when loading from cache
  x <- terra::writeRaster(x, filename = tempfile(fileext = ".tif"))

  return(list(x = x, class = "SpatRaster", unit = unit))
}
