#' toolProjectRegionalDelta
#'
#' Project amounts onto target data: for every year after the harmonization
#' year, the change of each region's amount total relative to the harmonization
#' year is distributed onto the cells of the region proportional to their share
#' weight of that year, on top of the per-cell amount anchor of the
#' harmonization year. Region totals are preserved, and amounts only grow where
#' the share weight is positive.
#'
#' @param xTarget target amounts, provides the per-cell anchor for the
#'   harmonization year
#' @param xInput input amounts, containing the harmonization year and the years
#'   after it
#' @param shareWeight weights defining each cell's share of its region's
#'   projected delta for the years after the harmonization year
#' @param harmonizationYear integer year the input changes are measured from
#' @return projected amounts for the years after the harmonization year,
#'   clamped at zero, in the item order of xTarget
#'
#' @author Pascal Sauer
toolProjectRegionalDelta <- function(xTarget, xInput, shareWeight, harmonizationYear) {
  yearsAfter <- getYears(xInput, as.integer = TRUE)
  yearsAfter <- yearsAfter[yearsAfter > harmonizationYear]
  regionAmount <- dimSums(xInput, 1.2)
  regionDelta <- regionAmount[, yearsAfter, ] - setYears(regionAmount[, harmonizationYear, ], NULL)
  share <- shareWeight / dimSums(shareWeight, 1.2)
  share[is.nan(share)] <- 0
  # the anchor comes first so that the result follows its item order, the
  # magclass arithmetic aligns the delta with it by name
  projected <- setYears(xTarget[, harmonizationYear, ], NULL) + regionDelta * share
  projected[projected < 0] <- 0
  stopifnot(identical(getItems(projected, 1), getItems(xTarget, 1)),
            identical(getItems(projected, 3), getItems(xTarget, 3)))
  return(projected)
}
