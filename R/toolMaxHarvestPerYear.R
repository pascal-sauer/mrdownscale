#' toolMaxHarvestPerYear
#'
#' Calculate the maximum possible harvest area per year
#' (dividing by timestep length) based on the given land data.
#'
#' @param land magpie object with at least the following categories:
#' c("primf", "secdf", "pltns", "primn", "secdn")
#' @param split if TRUE: split secdf to secyf and secdmf,
#' rename secdf to secnf, and append "_wood_harvest_area" to names
#' @return magpie object with the maximum possible yearly wood harvest area
#'
#' @author Pascal Sauer
toolMaxHarvestPerYear <- function(land, split = TRUE) {
  timestepLength <- new.magpie(years = getYears(land)[-1],
                               fill = diff(getYears(land, as.integer = TRUE)))
  stopifnot(timestepLength > 0)

  land <- land[, , c("primf", "secdf", "pltns", "primn", "secdn")]
  if (split) {
    getItems(land, 3) <- sub("secdf", "secyf", getItems(land, 3))
    getItems(land, 3) <- sub("secdn", "secnf", getItems(land, 3))
    getItems(land, 3) <- paste0(getItems(land, 3), "_wood_harvest_area")
    land <- add_columns(land, "secmf_wood_harvest_area")
    land[, , "secmf_wood_harvest_area"] <- land[, , "secyf_wood_harvest_area"]
  }
  maxHarvest <- setYears(land[, -nyears(land), ], getYears(land)[-1])
  maxHarvestPerYear <- maxHarvest / timestepLength
  return(maxHarvestPerYear)
}
