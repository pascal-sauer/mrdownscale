#' readLUH3
#'
#' Read LUH3 data. For the states subtype, the secma and secmb categories are removed.
#' For the management subtype, only the categories cpbf1, cpbf2, rndwd, fulwd, fertl and irrig are read.
#' For the transitions subtype, only the wood harvest categories bioh and harv are read. To
#' match magpie semantics years are shifted by 1 when reading transitions.
#' The LUH3 nc files have day-based time, which is converted to years.
#'
#' @param subtype one of states, management, transitions, cellArea
#' @param subset which years to read
#' @return data read from LUH3 historic nc files as SpatRaster
#'
#' @author Pascal Sauer
readLUH3 <- function(subtype, subset) {
  if (subtype == "cellArea") {
    cellArea <- terra::rast("multiple-static_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn.nc", "carea")
    return(list(x = cellArea, class = "SpatRaster", cache = FALSE, unit = "km2"))
  }

  years <- subset

  readLayers <- function(nc, variables, years) {
    yearIndizes <- years - 849 # LUH data starts in 850
    return(terra::rast(nc, lyrs = paste0(rep(variables, each = length(yearIndizes)), "_", yearIndizes)))
  }

  if (subtype == "states") {
    # all except secmb secma
    variables <- c("primf", "primn", "secdf", "secdn", "urban", "c3ann", "c4ann",
                   "c3per", "c4per", "c3nfx", "pastr", "range")
    x <- readLayers("multiple-states_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2024.nc",
                    variables, years)
    stopifnot(max(terra::minmax(x, compute = TRUE)) <= 1.0001,
              all(terra::units(x) == "1"))
    unit <- "1"
  } else if (subtype == "management") {
    cropTypes <- c("c3ann", "c4ann", "c3per", "c4per", "c3nfx")
    variables <- c(paste0("cpbf1_", cropTypes),
                   paste0("fertl_", cropTypes),
                   paste0("irrig_", cropTypes),
                   "flood", "rndwd", "fulwd", "combf")
    # the following variables are present but all zeros, so no need to read those:
    # cpbf2_*, pltns_wdprd, pltns_bfuel
    # prtct_primf, prtct_primn, prtct_secdf, prtct_secdn, prtct_pltns

    x <- readLayers("multiple-management_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2024.nc",
                    variables, years)

    # combf is a share of wood harvest like rndwd and fulwd, but we can ignore it as long as it is 0 everywhere
    stopifnot(all(terra::minmax(x["combf"], compute = TRUE) == 0))
    x <- x[[grep("combf", names(x), invert = TRUE)]]

    stopifnot(setequal(terra::varnames(x),
                       c(paste0("cpbf1_", cropTypes),
                         paste0("fertl_", cropTypes),
                         paste0("irrig_", cropTypes),
                         "flood", "rndwd", "fulwd")))

    unit <- "1, except fertl: kg ha-1 yr-1"
  } else if (subtype == "transitions") {
    woodland <- c("primf", "primn", "secmf", "secyf", "secnf", "pltns")
    variables <- c(paste0(woodland, "_harv"), paste0(woodland, "_bioh"))
    x <- readLayers("multiple-transitions_input4MIPs_landState_CMIP_UofMD-landState-3-1-1_gn_0850-2023.nc",
                    variables, years - 1)

    unit <- "*_bioh: kg C yr-1, *_harv: 1"
  } else {
    stop("subtype must be states, management, transitions or cellArea")
  }

  terra::time(x, tstep = "years") <- as.integer(sub("-01-01$", "", terra::time(x)))
  if (subtype == "transitions") {
    # LUH uses from-semantics for transitions (value for 1994 describes what happens from 1994 to 1995)
    # by adding 1 to time we get to-semantics (value for 1994 describes what happens from 1993 to 1994)
    terra::time(x, tstep = "years") <- terra::time(x) + 1
  }
  names(x) <- paste0("y", terra::time(x), "..", sub("_[0-9]+$", "", names(x)))

  return(list(x = x,
              class = "SpatRaster",
              cache = FALSE,
              unit = unit))
}
