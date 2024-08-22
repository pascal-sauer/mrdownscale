#' toolWoodHarvestArea
#'
#' Check wood harvest area is not exceeding land area of the corresponding
#' type in the previous timestep divided by timestep length. Optionally fix
#' by simply reducing wood harvest area to at most match land area.
#'
#' @param x magpie object with exactly the following categories:
#' paste0(c("primf", "secyf", "secmf", "primn", "secnf"), "_wood_harvest_area")
#' @param land magpie object with at least the following categories:
#' c("primf", "secdf", "primn", "secdn")
#' @param fix logical, if TRUE reduce wood harvest area to at most match land area
#' @return potentially fixed magpie object with the same dimensions as x
#'
#' @author Pascal Sauer
toolWoodHarvestArea <- function(x, land, fix) {
  stopifnot(setequal(getItems(x, 1), getItems(land, 1)))
  woodHarvestAreaBefore <- x

  woodHarvestArea <- .aggregateWoodHarvest(x)

  stopifnot(identical(getYears(woodHarvestArea), getYears(land)),
            getItems(woodHarvestArea, 3) %in% getItems(land, 3))
  years <- getYears(land, as.integer = TRUE)
  timestepLengths <- new.magpie(years = years[-1], fill = diff(years))
  woodland <- setYears(land[, -nyears(land), getItems(woodHarvestArea, 3)], years[-1])

  if (fix) {
    woodHarvestArea[, -1, ] <- pmin(woodHarvestArea[, -1, ], woodland / timestepLengths)
    x <- toolAggregate(woodHarvestArea, map, weight = woodHarvestAreaBefore + 10^-10,
                       from = "land", to = "harvest", dim = 3)

    maxCorrection <- max(abs(woodHarvestAreaBefore - x))
    if (maxCorrection > 10^-10) {
      toolStatusMessage("note", paste("wood harvest area was reduced to not exceed land area",
                                      "(max correction: ", signif(maxCorrection, 3), " Mha yr-1)"),
                        level = 1)
    }
  }

  landOvershoot <- min(woodland - timestepLengths * collapseDim(woodHarvestArea)[, -1, ])
  toolExpectTrue(landOvershoot >= -10^-10, paste0("Wood harvest area is smaller than land ",
                                                  "of the corresponding type ",
                                                  "(max overshoot: ", signif(-landOvershoot, 3), " Mha)"),
                 level = 1)
  return(x)
}

.aggregateWoodHarvest <- function(woodHarvest) {
  map <- as.data.frame(rbind(c("primf_wood_harvest_area", "primf"),
                              c("secyf_wood_harvest_area", "secdf"),
                              c("secmf_wood_harvest_area", "secdf"),
                              c("primn_wood_harvest_area", "primn"),
                              c("secnf_wood_harvest_area", "secdn")))
  colnames(map) <- c("harvest", "land")

  stopifnot(setequal(getItems(woodHarvest, 3), map$harvest))

  return(toolAggregate(woodHarvest, map, from = "harvest", to = "land", dim = 3))
}
