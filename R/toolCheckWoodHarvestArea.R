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
#' @param endOfHistory The last year considered part of the historical period,
#' will check and report consistency separately for history and after
#'
#' @author Pascal Sauer
toolCheckWoodHarvestArea <- function(harvest, land, endOfHistory) {
  stopifnot(identical(getItems(harvest, 1), getItems(land, 1)))
  if (nyears(harvest) < 2) {
    message("less than 2 years provided, not checking wood harvest area")
    return()
  }
  harvest <- toolAggregateWoodHarvest(harvest)
  stopifnot(getItems(harvest, 3) %in% getItems(land, 3),
            identical(getYears(harvest), getYears(land)))

  maxHarvestPerYear <- toolMaxHarvestPerYear(land, split = FALSE)
  excessHarvestPerYear <- harvest[, -1, ] - maxHarvestPerYear
  stopifnot(identical(getItems(excessHarvestPerYear, 2), getItems(maxHarvestPerYear, 2)))

  checkArea <- function(x, notePrefix) {
    if (nyears(x) == 0) {
      return()
    }
    msg <- paste0(" (max yearly excess harvest: ", signif(max(x), 3), " Mha)")
    toolExpectTrue(max(x) <= 10^-10,
                   paste0(notePrefix, "wood harvest area is smaller than land ",
                          "of the corresponding type", if (max(x) > 10^-10) msg))
  }
  harvestYears <- getYears(excessHarvestPerYear, as.integer = TRUE)
  checkArea(excessHarvestPerYear[, harvestYears[harvestYears <= endOfHistory], ],
            paste0("In historical period (until ", endOfHistory, "), "))
  checkArea(excessHarvestPerYear[, harvestYears[harvestYears > endOfHistory], ],
            paste0("After historical period (after ", endOfHistory, "), "))

  prim <- c("primf", "primn")
  landYears <- getYears(land, as.integer = TRUE)
  timestepLength <- new.magpie(years = landYears[-1], fill = diff(landYears))
  maxPrim <- setYears(land[, -nyears(land), prim], landYears[-1]) - timestepLength * harvest[, -1, prim]
  primExcess <- land[, -1, prim] - maxPrim

  checkPrim <- function(x, notePrefix) {
    if (nyears(x) == 0) {
      return()
    }
    msg <- paste0(" (", signif(max(x), 3), "Mha more primf/primn than possible)")
    toolExpectTrue(max(x) <= 10^-10,
                   paste0(notePrefix, "primf and primn are shrinking by at least ",
                          "their respective wood harvest area",
                          if (max(x) > 10^-10) msg))
  }
  primYears <- getYears(primExcess, as.integer = TRUE)
  checkPrim(primExcess[, primYears[primYears <= endOfHistory], ],
            paste0("In historical period (until ", endOfHistory, "), "))
  checkPrim(primExcess[, primYears[primYears > endOfHistory], ],
            paste0("After historical period (after ", endOfHistory, "), "))
}

toolWoodHarvestMapping <- function() {
  map <- as.data.frame(rbind(c("wood_harvest_area.primf", "primf"),
                             c("wood_harvest_area.secyf", "secdf"),
                             c("wood_harvest_area.secmf", "secdf"),
                             c("wood_harvest_area.pltns", "pltns"),
                             c("wood_harvest_area.primn", "primn"),
                             c("wood_harvest_area.secnf", "secdn")))
  colnames(map) <- c("harvest", "land")
  return(map)
}

toolAggregateWoodHarvest <- function(woodHarvest) {
  map <- toolWoodHarvestMapping()

  stopifnot(setequal(getItems(woodHarvest, 3), map$harvest))

  return(toolAggregate(woodHarvest, map, from = "harvest", to = "land", dim = 3))
}
