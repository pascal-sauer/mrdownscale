#' calcLandTarget
#'
#' Prepare the high resolution target land use dataset for
#' harmonization and downscaling, checking data for consistency before returning.
#'
#' @param target name of the target dataset, one of luh2, luh2mod, luh3
#' luh2mod/luh3 will split secdf into pltns and secdf
#' @param endOfHistory years later than this are not returned
#' @return land target data
#' @author Pascal Sauer
calcLandTarget <- function(target, endOfHistory) {
  x <- calcOutput("LandTargetComplete", target = target, aggregate = FALSE, supplementary = TRUE)
  x$x <- x$x[[terra::time(x$x) <= endOfHistory]]
  out <- x[c("x", "class", "unit", "description")]
  out$cache <- FALSE
  return(out)
}

calcLandTargetComplete <- function(target) {
  if (target %in% c("luh2", "luh2mod", "luh3")) {
    cropTypes <- c("c3ann", "c3nfx", "c3per", "c4ann", "c4per")

    if (target == "luh3") {
      states <- readSource("LUH3", subtype = "states", subset = 1995:2024)
    } else {
      states <- readSource("LUH2v2h", subtype = "states")
    }
    nonCropStates <- states[[grep(paste(cropTypes, collapse = "|"), names(states),
                                  value = TRUE, invert = TRUE)]]
    states <- spatRasterToDataset(states)

    if (target == "luh3") {
      man <- readSource("LUH3", subtype = "management", subset = 1995:2024, convert = FALSE)
      man <- man["irrig|cpbf1"]
      man <- spatRasterToDataset(man)
      cpbf1Category <- "cpbf1"
    } else {
      man <- readSource("LUH2v2h", subtype = "management", convert = FALSE)
      man <- spatRasterToDataset(man)
      man <- man[c(paste0("crpbf_", cropTypes),
                   paste0("irrig_", cropTypes))]
      cpbf1Category <- "crpbf"
    }
    stopifnot(all.equal(terra::time(states[1]), terra::time(man[1])))

    out <- do.call(c, c(lapply(cropTypes,
                               function(cropType) toolSplitIrrigationRainfed(man, cpbf1Category, cropType, states)),
                        nonCropStates))
    terra::time(out, tstep = "years") <- as.integer(substr(names(out), 2, 5))

    if (target %in% c("luh2mod", "luh3")) {
      out <- toolSplitSecdf(out)
    }
    expectedCategories <- unique(toolLandCategoriesMapping(input = "magpie", target = target)$dataOutput)
  } else if (target == "landuseinit") {
    out <- readSource("LanduseInit")
    getItems(out, 3) <- sub("primforest", "primf", getItems(out, 3))
    getItems(out, 3) <- sub("secdforest", "secdf", getItems(out, 3))

    out <- toolScaleConstantArea(out)
    out <- toolReplaceExpansion(out, "primf", "secdf")
    out <- as.SpatRaster(out)

    expectedCategories <- c("crop", "past", "forestry", "primf", "secdf", "urban",  "other")
  } else if (target == "landuseinitchina") {
    out <- readSource("LanduseInit")
    chinaCrops <- readSource("ChinaCrops")

    out <- out[getItems(chinaCrops, 1), , ]
    stopifnot(identical(getYears(out), getYears(chinaCrops)))

    otherCrop <- out[, , "crop"] - dimSums(chinaCrops, dim = 3)
    stopifnot(otherCrop >= 0)
    otherCrop <- magclass::setNames(otherCrop, "other_crop")

    out <- mbind(out[, , "crop", invert = TRUE], chinaCrops, otherCrop)

    getItems(out, 3) <- sub("primforest", "primf", getItems(out, 3))
    getItems(out, 3) <- sub("secdforest", "secdf", getItems(out, 3))

    out <- toolScaleConstantArea(out)
    out <- toolReplaceExpansion(out, "primf", "secdf")
    out <- as.SpatRaster(out)

    expectedCategories <- c("past", "forestry", "primforest", "secdforest", "urban", "other",
                            "rice_pro", "tece", "maiz", "other_crop")
  } else {
    stop("Unsupported output type \"", target, "\"")
  }

  # checks
  toolExpectTrue(terra::crs(out) != "", "Data contains CRS information")
  toolExpectTrue(setequal(sub("y[0-9]+\\.\\.", "", names(out)), expectedCategories),
                 "Land target categories match the corresponding mapping")
  toolExpectTrue(min(terra::minmax(out)) >= 0, "All values are >= 0")

  years <- sort(unique(terra::time(out)))
  primfn <- out["primf|primn"]
  primfnTime <- terra::time(primfn)
  primfnDiff <- primfn[[primfnTime %in% years[-1]]] - primfn[[primfnTime %in% years[-length(years)]]]
  toolExpectTrue(max(terra::minmax(primfnDiff)) <= 0,
                 "primary land is never expanding", falseStatus = "warn")

  return(list(x = out,
              class = "SpatRaster",
              unit = "Mha",
              description = "Land target data for data harmonization"))
}

toolSplitIrrigationRainfed <- function(man, cpbf1Category, cropType, states) {
  cpbf1 <- man[paste0(cpbf1Category, "_", cropType)]
  irrig <- man[paste0("irrig_", cropType)]
  stopifnot(all.equal(terra::time(cpbf1), terra::time(states[cropType])),
            all.equal(terra::time(irrig), terra::time(states[cropType])))

  # cpbf1 is 1st gen biofuel share of <cropType>
  # -> get <cropType>_biofuel_1st_gen area in Mha by multiplying with <cropType>
  biofuel1stGen <- cpbf1 * states[cropType]

  irrigatedBiofuel1stGen <- irrig * biofuel1stGen
  rainfedBiofuel1stGen <- biofuel1stGen - irrigatedBiofuel1stGen
  names(irrigatedBiofuel1stGen) <- sub("\\.\\..+$", paste0("..", cropType, "_irrigated_biofuel_1st_gen"),
                                       names(irrigatedBiofuel1stGen))
  names(rainfedBiofuel1stGen) <- sub("\\.\\..+$", paste0("..", cropType, "_rainfed_biofuel_1st_gen"),
                                     names(rainfedBiofuel1stGen))

  nonBiofuel <- states[cropType] - biofuel1stGen

  irrigatedBiofuel2ndGen <- NULL
  rainfedBiofuel2ndGen <- NULL
  if (cropType %in% c("c3per", "c4per")) {
    # 2nd gen biofuel is not part of LUH2v2h, but we need it for the harmonization, so fill with zeros
    # 2nd gen biofuel is part of LUH3, but it is all zeros (even over ocean), so instead of reading it...
    irrigatedBiofuel2ndGen <- 0 * states[cropType]
    rainfedBiofuel2ndGen <- 0 * states[cropType]
    names(irrigatedBiofuel2ndGen) <- paste0(names(irrigatedBiofuel2ndGen), "_irrigated_biofuel_2nd_gen")
    names(rainfedBiofuel2ndGen) <- paste0(names(rainfedBiofuel2ndGen), "_rainfed_biofuel_2nd_gen")
    nonBiofuel <- nonBiofuel - irrigatedBiofuel2ndGen - rainfedBiofuel2ndGen
  }

  irrigatedNonBiofuel <- irrig * nonBiofuel
  rainfedNonBiofuel <- nonBiofuel - irrigatedNonBiofuel
  names(irrigatedNonBiofuel) <- sub("\\.\\..+$", paste0("..", cropType, "_irrigated"), names(irrigatedNonBiofuel))
  names(rainfedNonBiofuel) <- sub("\\.\\..+$", paste0("..", cropType, "_rainfed"), names(rainfedNonBiofuel))

  out <- c(irrigatedNonBiofuel, rainfedNonBiofuel, irrigatedBiofuel1stGen,
           rainfedBiofuel1stGen, irrigatedBiofuel2ndGen, rainfedBiofuel2ndGen)
  # need to write raster to disk to avoid memory issues
  out <- terra::writeRaster(out, tempfile(fileext = ".tif"))
  return(out)
}

toolSplitSecdf <- function(x) {
  pltnsShare <- read.magpie(system.file("extdata/forestryShare.mz", package = "mrdownscale"))
  pltnsShare <- as.SpatRaster(pltnsShare)
  pltnsShare <- terra::crop(pltnsShare, x, extend = TRUE)

  bestFitYears <- vapply(unique(terra::time(x)), function(y) {
    # find the closest smaller available year, e.g. for 2019 use 2015
    return(Find(function(a) a <= y, terra::time(pltnsShare), right = TRUE))
  }, integer(1))
  pltns <- x["secdf"] * pltnsShare[[match(bestFitYears, terra::time(pltnsShare))]]

  names(pltns) <- sub("secdf", "pltns", names(pltns))
  stopifnot(terra::nlyr(x["secdf"]) == terra::nlyr(pltns))
  secdf <- x["secdf"] - pltns
  x <- c(x[[!grepl("secdf", names(x))]], pltns, secdf)

  # cannot cache SpatRaster with both in-memory and on-disk/file sources,
  # so write `x` to a tif file to get SpatRaster with a single source (the tif file)
  x <- terra::writeRaster(x, filename = tempfile(fileext = ".tif"))
  return(x)
}
