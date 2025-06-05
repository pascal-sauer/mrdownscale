#' calcNonlandInput
#'
#' Prepare the nonland input data for category mapping, checking data for consistency before returning.
#'
#' All "Land" functions deal with area data, as opposed to "Nonland" functions which deal with non-area
#' data such as the amount of applied fertilizer. These are treated differently, because for area
#' data other constraints apply, e.g. the total area must be constant over time.
#'
#' @param input name of an input dataset, currently only "magpie"
#' @return nonland input data
#' @author Pascal Sauer
calcNonlandInput <- function(input) {
  if (input == "magpie") {
    woodHarvestWeight <- readSource("MagpieFulldataGdx", subtype = "woodHarvestWeight")
    # convert from Pg DM yr-1 to kg C yr-1
    woodHarvestWeight <- woodHarvestWeight * 10^9 * 0.5

    woodHarvestWeightSource <- dimSums(woodHarvestWeight, dim = "woodType")
    woodHarvestWeightSource <- add_dimension(woodHarvestWeightSource, dim = 3.1,
                                             add = "category", "wood_harvest_weight")
    getSets(woodHarvestWeightSource)[["d3.2"]] <- "data"

    woodHarvestWeightType <- dimSums(woodHarvestWeight, dim = "source")
    woodHarvestWeightType <- add_dimension(woodHarvestWeightType, dim = 3.1,
                                           add = "category", "wood_harvest_weight_type")
    getSets(woodHarvestWeightType)[["d3.2"]] <- "data"

    woodHarvestArea <- readSource("MagpieFulldataGdx", subtype = "woodHarvestArea")
    woodHarvestArea <- dimSums(woodHarvestArea, dim = "ageClass")
    woodHarvestArea <- add_dimension(woodHarvestArea, dim = 3.1,
                                     add = "category", "wood_harvest_area")
    getSets(woodHarvestArea)[["d3.2"]] <- "data"

    # check wood harvest area * time step length (as unit is Mha yr-1) <=
    # land of the correponding type (in the previous timestep)
    land <- calcOutput("LandInput", input = input, aggregate = FALSE)
    stopifnot(identical(getYears(woodHarvestArea), getYears(land)))
    years <- getYears(land, as.integer = TRUE)
    timestepLengths <- new.magpie(years = years[-1], fill = diff(years))
    woodland <- setYears(land[, -nyears(land), getItems(woodHarvestArea, 3.2)], years[-1])
    toolExpectTrue(min(woodland - timestepLengths * collapseDim(woodHarvestArea)[, -1, ]) >= -10^-10,
                   "Wood harvest area is smaller than land of the corresponding type")

    fertilizer <- readSource("MagpieFulldataGdx", subtype = "fertilizer")
    fertilizer <- add_dimension(fertilizer, dim = 3.1,
                                add = "category", "fertilizer")
    getSets(fertilizer)[["d3.2"]] <- "data"
    # clusters without crop area are NA, replace with 0
    fertilizer[is.na(fertilizer)] <- 0
    # there are some negative values very close to zero, replace with 0
    stopifnot(min(fertilizer) >= -10^-10)
    fertilizer[fertilizer < 0] <- 0
    # convert from Tg yr-1 to kg yr-1
    fertilizer <- fertilizer * 10^9

    out <- mbind(woodHarvestWeightSource, woodHarvestWeightType, woodHarvestArea, fertilizer)
  } else {
    stop("Unsupported input dataset \"", input, "\"")
  }

  # check data for consistency
  toolExpectTrue(identical(unname(getSets(out)), c("region", "id", "year", "category", "data")),
                 "Dimensions are named correctly")
  toolExpectTrue(all(out >= 0), "All values are >= 0")

  return(list(x = out,
              isocountries = FALSE,
              unit = "harvest_weight: kg C yr-1; harvest_area: Mha yr-1; fertilizer: kg yr-1",
              min = 0,
              description = "Nonland input data for data harmonization and downscaling pipeline"))
}
