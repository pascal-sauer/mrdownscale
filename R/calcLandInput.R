#' calcLandInput
#'
#' Prepare the land input data for the category mapping, checking data for consistency before returning.
#' All "Land" functions deal with area data, as opposed to "Nonland" functions which deal with non-area
#' data such as the amount of applied fertilizer. These are treated differently, because for area
#' data other constraints apply, e.g. the total area must be constant over time.
#'
#' input = "magpie": includes the land use categories
#' past (pasture, including rangeland), forestry (managed forest plantations),
#' primforest, secdforest, urban, other (other land) and many
#' specific crop types. Furthermore, 1st gen biofuel is
#' added and filled with zeros. 1st gen biofuel is only modeled implicitly in
#' magpie via demand, and because of trade it is unclear on what area 1st gen
#' biofuel is grown, also 1st gen biofuel is quickly phased out in magpie, so
#' we fill biofuel_1st_gen with zeros and rely on the harmonization to produce
#' a plausible 1st gen biofuel time series.
#' input = "witch": includes a subset of LUH land use categories:
#' primf, primn, secdn, pastr, c4per, pltns, secdf, c3ann_irrigated, c3ann_rainfed
#' These are given as shares. A "rest" category is added so shares sum up to 1.
#'
#' @param input name of an input dataset, options: "magpie", "witch"
#' @return land input data
#' @author Jan Philipp Dietrich, Pascal Sauer
calcLandInput <- function(input) { # before adding args, consider: many functions @inheritParams from this function
  if (input == "magpie") {
    land <- readSource("MagpieFulldataGdx", subtype = "land")
    crop <- readSource("MagpieFulldataGdx", subtype = "crop")
    getItems(crop, dim = 3.1, full = TRUE) <- sub("\\.", "_", getItems(crop, dim = 3, full = TRUE))
    getItems(crop, dim = 3.2) <- NULL

    toolExpectLessDiff(land[, , "crop_area"],
                       dimSums(crop, dim = 3),
                       10^-5, "sum over all crops equals crop_area")
    land[, , "crop_area"] <- dimSums(crop, dim = 3) # ensure they are actually identical to avoid dividing by zero

    # in case we have no crop_area, but fallow and/or treecover: assign to bio energy trees
    fallowTreecover <- dimSums(land[, , c("crop_fallow", "crop_treecover")], 3)
    zeroCropArea <- collapseDim(land[, , "crop_area"]) == 0
    # keep fallow plus treecover only where crop_area is zero
    fallowTreecover <- fallowTreecover * ifelse(zeroCropArea, 1, 0)
    land[, , c("crop_fallow", "crop_treecover")] <- ifelse(zeroCropArea, 0,
                                                           land[, , c("crop_fallow", "crop_treecover")])
    crop[, , "betr_rainfed"] <- crop[, , "betr_rainfed"] + fallowTreecover
    land[, , "crop_area"] <- land[, , "crop_area"] + fallowTreecover

    # scale crop to take up the whole area of crop_area + crop_fallow + crop_treecover
    # need this to report for ScenarioMIP/LUH-format, might not want this for other applications
    totalCrop <- dimSums(land[, , c("crop_area", "crop_fallow", "crop_treecover")], 3)
    scalingFactors <- totalCrop / dimSums(crop, dim = 3)
    scalingFactors[is.nan(scalingFactors)] <- 1
    stopifnot(1 <= scalingFactors, scalingFactors < Inf)
    crop <- crop * scalingFactors

    toolExpectLessDiff(dimSums(crop, dim = 3), totalCrop, 10^-5,
                       "after scaling, sum over all crops equals crop_area + crop_fallow + crop_treecover")

    out <- mbind(land[, , c("crop_area", "crop_fallow", "crop_treecover"), invert = TRUE], crop)

    # see note in the documentation of this function
    out <- add_columns(out, "biofuel_1st_gen", fill = 0)

    expectedCategories <- unique(toolGetMapping("referenceMappings/magpie.csv", where = "mrdownscale")$data)
    primf <- "primforest"
  } else if (input == "witch") {
    x <- readSource("WITCH")
    stopifnot(length(unique(x$model)) == 1,
              length(unique(x$scenario)) == 1,
              "world" %in% x$region)
    x <- x[x$region != "world", ]

    stopifnot(x[x$variable == "Land Cover", ]$units == "million ha")
    regionAreaMha <- x[x$variable == "Land Cover", c("region", "value")]
    stopifnot(length(unique(regionAreaMha$region)) == nrow(unique(regionAreaMha)))
    regionAreaMha <- unique(regionAreaMha)
    regionAreaMha <- collapseDim(as.magpie(regionAreaMha, spatial = "region"))

    x <- x[x$variable %in% c("primf", "secdf", "pltns", "primn", "secdn", "pastr",
                             "c3ann", "irrig_c3ann", "c4per"), ]

    # check irrig_c3ann is fraction of crop area, need to handle that later
    stopifnot(identical(unique(x[x$variable == "irrig_c3ann", ]$units),
                        "fraction of crop area"))
    x[x$variable == "irrig_c3ann", ]$units <- "fraction of grid cell"
    stopifnot(identical(unique(x$units), "fraction of grid cell"))

    x <- x[, c("region", "year", "variable", "value")]
    out <- as.magpie(x, spatial = "region", temporal = "year")

    # add artificial region numbers/ids as these are expected later
    mapping <- readSource("WITCH", subtype = "resolutionMapping")
    out <- toolAggregate(out, unique(mapping[, c("witch17", "lowRes")]))
    names(dimnames(out))[1] <- "region.id"

    if (anyNA(out)) {
      toolStatusMessage("warn", "NAs detected, replacing with 0.")
      out[is.na(out)] <- 0
    }
    if (min(out) < 0) {
      toolStatusMessage("warn", "Negative values detected, replacing with 0.")
      # this leads to sum of shares > 1, scaling below won't be needed once this is fixed
      out[out < 0] <- 0
    }
    names(dimnames(out))[3] <- "data"

    out <- add_columns(out, c("c3ann_irrigated", "c3ann_rainfed"))

    # convert from fraction of crop area to fraction of grid cell
    stopifnot(out[, , "irrig_c3ann"] <= 1)
    out[, , "c3ann_irrigated"] <- out[, , "irrig_c3ann"] * out[, , "c3ann"]

    out[, , "c3ann_rainfed"] <- (1 - out[, , "irrig_c3ann"]) * out[, , "c3ann"]
    out <- out[, , c("c3ann", "irrig_c3ann"), invert = TRUE]

    shareTotal <- dimSums(out, 3)
    if (max(shareTotal - 1) > 10^-10) {
      toolStatusMessage("warn", "Scaling land so that land share sum is equal to 1.")
      shareTotal[shareTotal < 1] <- 1
      out <- out / shareTotal
    }

    rest <- magclass::setNames(1 - dimSums(out, 3), "rest")
    rest[rest < 0] <- 0
    out <- mbind(out, rest)

    toolExpectLessDiff(dimSums(out, 3), 1, 10^-10, "land shares sum up to 1")

    # convert from fraction of grid cell to Mha
    out <- out * regionAreaMha

    expectedCategories <- c("primf", "secdf", "pltns", "primn", "secdn", "pastr",
                            "c3ann_irrigated", "c3ann_rainfed", "c4per", "rest")
    primf <- "primf"
  } else if (input == "coffee") {
    x <- readSource("COFFEE")

    stopifnot(length(unique(x$Model)) == 1,
              length(unique(x$Scenario)) == 1,
              "World" %in% x$Region,
              x$Unit == "Mha")

    if (min(x$Value) < 0) {
      toolStatusMessage("warn", "Negative values detected, replacing with 0.")
      negative <- x[x$Value < 0, ]
      print(utils::head(negative[order(negative$Value), ], n = 5))
      x$Value[x$Value < 0] <- 0
    }

    x <- x[x$Region != "World", ]

    x <- x[, c("Region", "Year", "Variable", "Value")]
    out <- as.magpie(x, spatial = "Region", temporal = "Year")

    # add artificial region numbers/ids as these are expected later
    mapping <- readSource("COFFEE", subtype = "regionMapping", convert = FALSE)
    out <- toolAggregate(out, unique(mapping[, c("region", "lowRes")]))
    names(dimnames(out)) <- c("region.id", "year", "data")

    refmap <- toolGetMapping("referenceMappings/coffee.csv", where = "mrdownscale")
    vars <- setdiff(unique(refmap$data), "rest")
    if (any(out[, , "Land_Cover"] - dimSums(out[, , vars], 3) < 0)) {
      warning("Land Cover is smaller than sum of other variables, ",
              "so multiply by 1.1 until this is fixed in the input data")
      rest <- magclass::setNames(1.1 * out[, , "Land_Cover"] - dimSums(out[, , vars], 3), "rest")
    } else {
      rest <- magclass::setNames(out[, , "Land_Cover"] - dimSums(out[, , vars], 3), "rest")
    }
    stopifnot(rest >= 0)
    out <- mbind(out, rest)

    out <- out[, , unique(refmap$data)]

    expectedCategories <- unique(refmap$data)
    primf <- "Land_Cover_Forest_Primary"
  } else {
    stop("Unsupported input type \"", input, "\"")
  }

  # check data for consistency
  toolExpectTrue(identical(unname(getSets(out)), c("region", "id", "year", "data")),
                 "Dimensions are named correctly")
  toolExpectTrue(setequal(getItems(out, dim = 3), expectedCategories),
                 "Land input categories match expactation")
  toolExpectTrue(all(out >= 0), "All values are >= 0")
  outSum <- dimSums(out, dim = 3)
  toolExpectLessDiff(outSum, outSum[, 1, ], 10^-4, "Total area is constant over time")
  toolExpectTrue(all(out[, -1, primf] <= setYears(out[, -nyears(out), primf], getYears(out[, -1, ]))),
                 "primary forest is never expanding", falseStatus = "warn")

  return(list(x = out,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Land input data for data harmonization and downscaling pipeline",
              clean_magpie = FALSE)) # preserve region ids
}
