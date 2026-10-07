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
#' compensated by the other categories of the group. For crops, negative values
#' are first compensated by the corresponding rainfed/irrigated twin of the
#' same crop, which keeps that pair's total area, including biofuel types like
#' c3ann_rainfed_biofuel_1st_gen. Crop categories without their twin are
#' compensated by the whole group. Groups with a negative total area are set to
#' zero. Afterwards, if any negatives remain, all categories except urban are
#' scaled down to keep the total area constant.
#' This group compensation and constant total area applies only to land use
#' data, which is identified by its category taxonomy. For nonland data, which
#' is not area conservative, negative values are simply clamped to zero.
#'
#' @param xInput input data as magpie object
#' @param xTarget target data as magpie object
#' @param harmonizationPeriod Two identical integer values, the year the
#'   absolute changes of the input data are applied to, must be present in both
#'   input and target data
#' @return harmonized data set as magpie object with data from target for years
#'   up to and including the harmonization year and absolute changes from input
#'   relative to the harmonization year afterwards
#' @author Pascal Sauer
toolHarmonizeAbsoluteChanges <- function(xInput, xTarget, harmonizationPeriod) {
  hp <- harmonizationPeriod
  stopifnot(length(hp) == 2,
            round(hp) == hp,
            identical(hp[1], hp[2]))
  hy <- hp[1]

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)

  stopifnot(!anyNA(xInput),
            !anyNA(xTarget),
            setequal(getItems(xInput, 1), getItems(xTarget, 1)),
            setequal(getItems(xInput, 3), getItems(xTarget, 3)),
            hy %in% inputYears,
            hy %in% targetYears)
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  # apply absolute changes of input data to target data of the harmonization year
  # negative numbers are possible if input data loses more area of a
  # category than the target data has in the harmonization year
  raw <- setYears(xTarget[, hy, ], NULL) + (xInput[, inputYears > hy, ] - setYears(xInput[, hy, ], NULL))
  changed <- raw

  croplandPattern <- "_rainfed|_irrigated"
  groups <- list(
    "forest and other land" = intersect(c("pltns", "primf", "secdf", "primn", "secdn"), getItems(changed, 3)),
    cropland = grep(croplandPattern, getItems(changed, 3), value = TRUE),
    "pasture and rangeland" = intersect(c("pastr", "range"), getItems(changed, 3))
  )
  # like the fade harmonizers land use data is identified by its category
  # taxonomy, but also partial overlap counts as land so that data with
  # invalid categories is rejected instead of treated as nonland
  isLand <- any(getItems(changed, 3) %in% c(unlist(groups), "urban"))

  if (isLand) {
    # total area of each cell, constant over time
    targetArea <- dimSums(setYears(xTarget[, hy, ], NULL), 3)
    stopifnot(all(abs(dimSums(xTarget, 3) - targetArea) < 10^-5))
    stopifnot(all(abs(dimSums(changed, 3) - targetArea) < 10^-5))

    stopifnot(setequal(c(unlist(groups), "urban"), getItems(changed, 3)))
    # for crops, negative values are first compensated by scaling
    # the corresponding rainfed/irrigated twin of the same crop.
    # Categories without twin (incl. all non-crops) are compensated by scaling the whole group.
    for (group in groups) {
      groupArea <- changed[, , group]
      groupTarget <- pmax(dimSums(groupArea, 3), 0)

      cropItems <- grep(croplandPattern, group, value = TRUE)
      pairKey <- sub(croplandPattern, "", cropItems)
      pairs <- split(cropItems, pairKey)
      for (pair in pairs[lengths(pairs) == 2]) {
        pairArea <- groupArea[, , pair]
        groupArea[, , pair] <- toolHandleNegatives(pairArea, targetArea = pmax(dimSums(pairArea, 3), 0))
      }

      changed[, , group] <- toolHandleNegatives(groupArea, targetArea = groupTarget)
    }

    changed[, , "urban"] <- pmax(changed[, , "urban"], 0)

    nonUrbanTarget <- targetArea - dimSums(changed[, , "urban"], 3)
    stopifnot(all(nonUrbanTarget >= -10^-5))
    nonUrbanTarget[nonUrbanTarget < 0] <- 0
    nonUrban <- setdiff(getItems(changed, 3), "urban")
    changed[, , nonUrban] <- toolHandleNegatives(changed[, , nonUrban], targetArea = nonUrbanTarget)

    out <- mbind(xTarget[, targetYears <= hy, ], changed)

    # prim expansion is expected after harmonization due to prim differences between input and target dataset
    out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100)
    out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

    toolReportAreaDeviation(raw, out, groups = groups)

    stopifnot(all(abs(dimSums(out, 3) - targetArea) < 10^-5))
  } else {
    # nonland data is not area conservative, apply the raw absolute changes to
    # the target data and just clamp negative values to zero
    out <- mbind(xTarget[, targetYears <= hy, ], raw)
    out[out < 0] <- 0
  }

  return(out)
}
