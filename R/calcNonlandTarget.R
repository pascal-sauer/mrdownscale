#' calcNonlandTarget
#'
#' Prepare the high resolution nonland target dataset for
#' harmonization and downscaling, checking data for consistency before returning.
#'
#' @param target name of a target dataset, see \code{\link{calcLandTarget}}
#' for available target datasets
#' @param endOfHistory years later than this are not returned
#' @return nonland target data
#' @author Pascal Sauer
calcNonlandTarget <- function(target, endOfHistory) {
  x <- calcOutput("NonlandTargetComplete", target = target, aggregate = FALSE, supplementary = TRUE)
  x$x <- x$x[[terra::time(x$x) <= endOfHistory]]
  out <- x[c("x", "class", "unit", "description")]
  out$cache <- FALSE
  return(out)
}

calcNonlandTargetComplete <- function(target) {
  if (target %in% c("luh2", "luh2mod", "luh3pltns")) {
    if (target %in% c("luh2", "luh2mod")) {
      cellAreaKm2 <- readSource("LUH2v2h", subtype = "cellArea", convert = FALSE)
      management <- readSource("LUH2v2h", subtype = "management", convert = FALSE)
      transitions <- readSource("LUH2v2h", subtype = "transitions", convert = FALSE)
    } else {
      cellAreaKm2 <- readSource("LUH3", subtype = "cellArea", convert = FALSE)
      management <- readSource("LUH3", subtype = "management", subset = 1995:2024, convert = FALSE)
      transitions <- readSource("LUH3", subtype = "transitions", subset = 1995:2024, convert = FALSE)
    }

    # convert from km2 to Mha
    cellAreaMha <- cellAreaKm2 / 10000

    fertilizer <- management["fertl"]
    names(fertilizer) <- paste0(sub("fertl_", "", names(fertilizer)), "_fertilizer")

    ### wood harvest area in Mha yr-1
    # LUH3 includes pltns_harv variable, but it's zero everywhere, so no need to read it
    pltnsWoodHarvestArea <- transitions["primf_harv"] * 0
    names(pltnsWoodHarvestArea) <- sub("primf", "pltns", names(pltnsWoodHarvestArea))
    terra::varnames(pltnsWoodHarvestArea) <- "pltns_harv"
    terra::longnames(pltnsWoodHarvestArea) <- "wood harvest area from plantation forest vegetation"

    # convert from shares to Mha yr-1
    woodHarvestArea <- c(transitions["primf_harv"] * cellAreaMha,
                         transitions["primn_harv"] * cellAreaMha,
                         transitions["secmf_harv"] * cellAreaMha,
                         transitions["secyf_harv"] * cellAreaMha,
                         transitions["secnf_harv"] * cellAreaMha,
                         pltnsWoodHarvestArea)
    names(woodHarvestArea) <- paste0(sub("_harv", "", names(woodHarvestArea)), "_wood_harvest_area")
    terra::units(woodHarvestArea) <- "Mha yr-1"

    ### wood harvest weight (bioh) in kg C yr-1
    woodHarvestWeight <- c(transitions["(primf|primn|secmf|secyf|secnf)_bioh"])

    # LUH3 includes pltns_bioh variable, but it's zero everywhere, so no need to read it
    pltnsBioh <- woodHarvestWeight["primf_bioh"] * 0
    names(pltnsBioh) <- sub("primf", "pltns", names(pltnsBioh))
    terra::varnames(pltnsBioh) <- "pltns_bioh"
    terra::longnames(pltnsBioh) <- "wood harvest biomass carbon from plantation forest vegetation"

    woodHarvestWeight <- c(woodHarvestWeight, pltnsBioh)

    minWoodHarvestWeight <- min(terra::minmax(woodHarvestWeight, compute = TRUE))
    if (minWoodHarvestWeight < 0) {
      # replace negative weight of wood harvest with 0
      toolStatusMessage("note", paste0("replacing negative wood harvest weight (bioh) with 0 (min: ",
                                       round(minWoodHarvestWeight, 3), " kg C yr-1)"))
      woodHarvestWeight <- terra::classify(woodHarvestWeight, cbind(-Inf, 0, 0))
    }
    terra::units(woodHarvestWeight) <- "kg C yr-1"

    ### wood harvest weight type (fuelwood/roundwood) in kg C yr-1
    years <- unique(terra::time(woodHarvestWeight))
    woodHarvestWeightType <- do.call(c, lapply(years, function(year) {
      total <- sum(woodHarvestWeight[[terra::time(woodHarvestWeight) == year]])
      roundwood <- total * management[[paste0("y", year, "..rndwd")]]
      names(roundwood) <- paste0("y", year, "..roundwood_harvest_weight_type")
      fuelwood <- total * management[[paste0("y", year, "..fulwd")]]
      names(fuelwood) <- paste0("y", year, "..fuelwood_harvest_weight_type")
      return(c(roundwood, fuelwood))
    }))
    terra::units(woodHarvestWeightType) <- "kg C yr-1"

    out <- c(woodHarvestArea, woodHarvestWeight, woodHarvestWeightType, fertilizer)
    terra::time(out, tstep = "years") <- as.integer(sub("^y([0-9]+).+", "\\1", names(out)))

    # cannot cache SpatRaster with both in-memory and on-disk/file sources,
    # so write `out` to a tif file to get SpatRaster with a single source (the tif file)
    out <- terra::writeRaster(out, filename = tempfile(fileext = ".tif"))

    toolExpectTrue(min(terra::minmax(out)) >= 0, "All values are >= 0")

    return(list(x = out,
                class = "SpatRaster",
                unit = "harvest_weight & bioh: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg ha-1 yr-1",
                description = "Nonland target data for data harmonization"))
  } else {
    stop("Unsupported output type \"", target, "\"")
  }
}
