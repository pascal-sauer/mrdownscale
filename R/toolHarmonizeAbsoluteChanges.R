#' toolHarmonizeAbsoluteChanges
#'
#' Tool function for creating a harmonized data set by applying the absolute
#' changes of the input data to the target data: up to the harmonization year
#' the target data is used, afterwards the difference of the input data to the
#' input data of the harmonization year is added to the target data of the
#' harmonization year.
#'
#' Negative values are handled within groups of related categories (forest &
#' other land; all cropland types; pasture & rangeland), so that they are
#' compensated by the other categories of the group. Groups with a negative
#' total area are set to zero. Afterwards, if any negatives remain, all
#' categories except urban are scaled down to keep the total area constant.
#'
#' @param xInput input data as magpie object
#' @param xTarget target data as magpie object
#' @param harmonizationPeriod Single integer value, the year the absolute
#'   changes of the input data are applied to, must be present in both input
#'   and target data
#' @return harmonized data set as magpie object with data from target for years
#'   up to and including the harmonization year and absolute changes from input
#'   relative to the harmonization year afterwards
#' @author Pascal Sauer
toolHarmonizeAbsoluteChanges <- function(xInput, xTarget, harmonizationPeriod) {
  hy <- harmonizationPeriod

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)

  stopifnot(length(hy) == 1,
            round(hy) == hy,
            !anyNA(xInput),
            !anyNA(xTarget),
            setequal(getItems(xInput, 1), getItems(xTarget, 1)),
            setequal(getItems(xInput, 3), getItems(xTarget, 3)),
            hy %in% inputYears,
            hy %in% targetYears)
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  # total area of each cell, constant over time
  targetArea <- dimSums(setYears(xTarget[, hy, ], NULL), 3)
  stopifnot(all(abs(dimSums(xTarget, 3) - targetArea) < 10^-5))

  # apply absolute changes of input data to target data of the harmonization year
  changed <- setYears(xTarget[, hy, ], NULL) + (xInput[, inputYears > hy, ] - setYears(xInput[, hy, ], NULL))
  stopifnot(all(abs(dimSums(changed, 3) - targetArea) < 10^-5))

  # absolute changes can become negative if the input data loses more area of a
  # category than the target data has in the harmonization year
  groups <- list(
    forestOther = intersect(c("pltns", "primf", "secdf", "primn", "secdn"), getItems(changed, 3)),
    cropland = grep("rainfed|irrigated", getItems(changed, 3), value = TRUE),
    pastureRangeland = intersect(c("pastr", "range"), getItems(changed, 3))
  )
  # set negatives to zero, then scale other variables from that group to achieve target area
  for (group in groups) {
    changed[, , group] <- toolHandleNegatives(changed[, , group], allowNegativeTarget = TRUE)
  }

  urban <- changed[, , "urban"]
  stopifnot(all(urban >= -10^-5))
  urban[urban < 0] <- 0
  changed[, , "urban"] <- urban
  nonUrbanTarget <- targetArea - dimSums(urban, 3)
  stopifnot(all(nonUrbanTarget >= 0))
  nonUrban <- setdiff(getItems(changed, 3), "urban")
  changed[, , nonUrban] <- toolHandleNegatives(changed[, , nonUrban], targetArea = nonUrbanTarget)

  out <- mbind(xTarget[, targetYears <= hy, ], changed)

  # prim expansion is expected after harmonization due to prim differences between input and target dataset
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

  stopifnot(all(abs(dimSums(out, 3) - targetArea) < 10^-5))

  return(out)
}
