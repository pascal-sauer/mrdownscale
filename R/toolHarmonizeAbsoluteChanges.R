#' toolHarmonizeAbsoluteChanges
#'
#' Tool function for creating a harmonized data set by applying the absolute
#' changes of the input data to the target data: up to the harmonization year
#' the target data is used, afterwards the difference of the input data to the
#' input data of the harmonization year is added to the target data of the
#' harmonization year.
#' If this makes forest categories (pltns, primf, secdf) negative, the missing
#' area is deducted proportionally from the other land categories (primn,
#' secdn). Remaining negative values are set to 0 and categories are scaled
#' down to keep the total area constant: first everything except
#' prim (primf, primn) and urban, then prim only if prim + urban exceed
#' total area. Urban is never scaled.
#'
#' @param xInput input data as magpie object
#' @param xTarget target data as magpie object
#' @param harmonizationPeriod Single integer value, the year the absolute
#' changes of the input data are applied to, must be present in both input and
#' target data
#' @return harmonized data set as magpie object with data from target for years
#' up to and including the harmonization year and absolute changes from input
#' relative to the harmonization year afterwards
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

  forest <- intersect(c("pltns", "primf", "secdf"), getItems(changed, 3))
  otherLand <- intersect(c("primn", "secdn"), getItems(changed, 3))

  # absolute changes can become negative if the input data loses more area of a
  # category than the target data has in the harmonization year
  # handle negative forest by reducing otherLand instead
  if (any(changed[, , forest] < 0)) {
    negatives <- as.vector(changed[, , forest])
    negatives <- negatives[negatives < 0]
    toolStatusMessage("note", paste0("absolute changes made forest categories negative: ",
                                     length(negatives), " of ", length(changed), " cluster-year-combinations (",
                                     round(100 * length(negatives) / length(changed), 1),
                                     "%), min = ", min(negatives),
                                     ", mean = ", mean(negatives),
                                     ", median = ", median(negatives),
                                     ", shortfall was deducted from ",
                                     paste(otherLand, collapse = "/")))
    # deduct missing area from other land categories instead of letting forest
    # categories become negative, proportional to their current shares
    shortfallForest <- -changed[, , forest] * (changed[, , forest] < 0)
    stopifnot(all(shortfallForest >= 0))
    totalShortfallForest <- dimSums(shortfallForest, 3)

    availOther <- changed[, , otherLand] * (changed[, , otherLand] > 0)
    stopifnot(all(availOther >= 0))
    totalAvailOther <- dimSums(availOther, 3)

    deducted <- pmin(totalShortfallForest, totalAvailOther)

    fact <- shortfallForest / (totalShortfallForest + (totalShortfallForest == 0))
    stopifnot(all(0 <= fact & fact <= 1))
    changed[, , forest] <- changed[, , forest] + deducted * fact

    fact <- availOther / (totalAvailOther + (totalAvailOther == 0))
    stopifnot(all(0 <= fact & fact <= 1))
    changed[, , otherLand] <- changed[, , otherLand] - deducted * fact
  }

  # set negative values to 0 and scale to match targetArea
  if (any(changed < 0)) {
    changed[changed < 0] <- 0

    stopifnot(all(dimSums(changed, 3) >= targetArea - 10^-5))

    prim <- intersect(c("primf", "primn"), getItems(changed, 3))
    rest <- setdiff(getItems(changed, 3), c(prim, "urban"))

    # scale items so that together with protectedSum they fit into the targetArea
    .reduceToFit <- function(items, protected) {
      itemSum <- dimSums(changed[, , items], 3)
      protectedSum <- dimSums(changed[, , protected], 3)
      factor <- (targetArea - protectedSum) / (itemSum + (itemSum == 0))
      factor[factor < 0] <- 0
      factor[factor > 1] <- 1
      return(changed[, , items] * factor)
    }

    changed[, , rest] <- .reduceToFit(rest, c(prim, "urban"))
    if (any(dimSums(changed[, , c(prim, "urban")], 3) > targetArea + 10^-5)) {
      toolStatusMessage("warn", paste0("prim + urban exceed total area after correcting ",
                                       "negative values, reducing prim categories"))
      changed[, , prim] <- .reduceToFit(prim, "urban")
      if (any(dimSums(changed[, , "urban"], 3) > targetArea + 10^-5)) {
        stop("urban area exceeds total area after correcting negative values, ",
             "urban would need to be scaled but that is not allowed")
      }
    }
  }

  out <- mbind(xTarget[, targetYears <= hy, ], changed)

  # prim expansion is expected after harmonization due to prim differences between input and target dataset
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

  stopifnot(all(abs(dimSums(out, 3) - targetArea) < 10^-5))

  return(out)
}
