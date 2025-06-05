#' toolCheckWoodHarvestArea
#'
#' Check wood harvest area is not exceeding land area of the corresponding
#' type divided by timestep length. Also, check that primf and primn
#' are reduced by at least as much as they were harvested.
#'
#' @param harvest magpie object with exactly the following categories:
#' paste0(c("primf", "secyf", "secmf", "pltns", "primn", "secnf"), "_wood_harvest_area")
#' @param land magpie object with at least the following categories:
#' c("primf", "secdf", "pltns", "primn", "secdn")
#' @param notePrefix character to prepend to the check's message
#'
#' @author Pascal Sauer
toolCheckWoodHarvestArea <- function(harvest, land, notePrefix = "") {
  stopifnot(setequal(getItems(harvest, 1), getItems(land, 1)))
  harvest <- toolAggregateWoodHarvest(harvest)
  stopifnot(getItems(harvest, 3) %in% getItems(land, 3),
            identical(getYears(harvest), getYears(land)))
  years <- getYears(land, as.integer = TRUE)

  maxHarvestPerYear <- toolMaxHarvestPerYear(land, split = FALSE)
  excessHarvestPerYear <- harvest[, -1, ] - maxHarvestPerYear
  stopifnot(identical(getItems(excessHarvestPerYear, 2), getItems(maxHarvestPerYear, 2)))

  msg <- paste0(" (max yearly excess harvest: ", signif(max(excessHarvestPerYear), 3), " Mha)")
  toolExpectTrue(max(excessHarvestPerYear) <= 10^-10,
                 paste0(notePrefix, "wood harvest area is smaller than land ",
                        "of the corresponding type", if (max(excessHarvestPerYear) > 10^-10) msg),
                 level = 1)

  prim <- c("primf", "primn")
  timestepLength <- new.magpie(years = years[-1], fill = diff(years))
  maxPrim <- setYears(land[, -nyears(land), prim], years[-1]) - timestepLength * harvest[, -1, prim]
  primExcess <- land[, -1, prim] - maxPrim

  msg <- paste0(" (", signif(max(primExcess), 3), "Mha more primf/primn than possible)")
  toolExpectTrue(max(primExcess) <= 10^-10,
                 paste0(notePrefix, "primf and primn are shrinking by at least ",
                        "their respective wood harvest area",
                        if (max(primExcess) > 10^-10) msg),
                 level = 1)
}

toolWoodHarvestMapping <- function() {
  map <- as.data.frame(rbind(c("primf_wood_harvest_area", "primf"),
                             c("secyf_wood_harvest_area", "secdf"),
                             c("secmf_wood_harvest_area", "secdf"),
                             c("pltns_wood_harvest_area", "pltns"),
                             c("primn_wood_harvest_area", "primn"),
                             c("secnf_wood_harvest_area", "secdn")))
  colnames(map) <- c("harvest", "land")
  return(map)
}

toolBiohMapping <- function() {
  map <- as.data.frame(rbind(c("primf_bioh", "primf"),
                             c("secyf_bioh", "secdf"),
                             c("secmf_bioh", "secdf"),
                             c("pltns_bioh", "pltns"),
                             c("primn_bioh", "primn"),
                             c("secnf_bioh", "secdn")))
  colnames(map) <- c("bioh", "land")
  return(map)
}

toolAggregateWoodHarvest <- function(woodHarvest) {
  map <- toolWoodHarvestMapping()

  stopifnot(setequal(getItems(woodHarvest, 3), map$harvest))

  return(toolAggregate(woodHarvest, map, from = "harvest", to = "land", dim = 3))
}

woodHarvestAreaCategories <- function() {
  return(paste0(c("primf", "secyf", "secmf", "pltns", "primn", "secnf"), "_wood_harvest_area"))
}
