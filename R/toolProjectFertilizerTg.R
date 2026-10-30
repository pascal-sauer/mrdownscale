#' toolProjectFertilizerTg
#'
#' Project fertilizer tonnages onto the harmonized land: for every year after
#' the harmonization year, the change of each region's tonnage total relative to
#' the harmonization year is distributed onto the clusters of the region
#' proportional to their harmonized crop group area of that year, on top of the
#' per-cluster target tonnage of the harmonization year. Region totals are
#' preserved, and tonnes only grow where the harmonized land has cropland of
#' the respective crop group.
#'
#' @param fertilizerTgTarget target fertilizer in Tg yr-1, provides the per-cluster
#'   tonnage anchor for the harmonization year
#' @param fertilizerTgInput input fertilizer in Tg yr-1, containing the
#'   harmonization year and the years after it
#' @param landMha harmonized land in Mha
#' @param harmonizationYear integer year the input changes are measured from
#' @return projected fertilizer in Tg yr-1 for the years after the
#'   harmonization year, clamped at zero, in the item order of
#'   fertilizerTgTarget
#'
#' @author Pascal Sauer
toolProjectFertilizerTg <- function(fertilizerTgTarget, fertilizerTgInput, landMha, harmonizationYear) {
  yearsAfter <- getYears(fertilizerTgInput, as.integer = TRUE)
  yearsAfter <- yearsAfter[yearsAfter > harmonizationYear]
  areaMha <- toolAggregateCropland(landMha[getItems(fertilizerTgTarget, 1), yearsAfter, ],
                                   getItems(fertilizerTgInput, 3.2), keepOthers = FALSE)
  return(toolProjectRegionalDelta(fertilizerTgTarget, fertilizerTgInput, areaMha, harmonizationYear))
}
