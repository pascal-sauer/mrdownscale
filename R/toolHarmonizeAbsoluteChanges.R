#' toolHarmonizeAbsoluteChanges
#'
#' Tool function for creating a harmonized data set by applying the absolute
#' changes of the input data to the target data. Before and in the
#' harmonization year the target data is used. Afterwards, for every year the
#' difference of the input data to the input data in the harmonization year is
#' added to the target data of the harmonization year, e.g. if the target data
#' has 10 Mha secdf in the harmonization year and the input data secdf grows by
#' 2 Mha from the harmonization year to 2025, the harmonized data set has 12
#' Mha secdf in 2025.
#' If this makes forest categories (pltns, primf, secdf) negative, the missing
#' area is deducted from the other land categories (primn, secdn) instead.
#' Negative values which remain afterwards (because other land was not
#' available in sufficient amount or non-forest categories became negative) are
#' set to 0 and the remaining categories are then scaled down, so that the total
#' area stays constant. In this process the other categories (everything except
#' primf, primn and urban) are scaled first and primf and primn only if prim and
#' urban together still exceed the total area. Urban is never scaled; if urban
#' alone exceeds the total area, an error is raised.
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
            setequal(getItems(xInput, dim = 1), getItems(xTarget, dim = 1)),
            setequal(getItems(xInput, dim = 3), getItems(xTarget, dim = 3)))
  if (!hy %in% targetYears) {
    stop("harmonizationPeriod ", hy, " is not a year in the target data")
  }
  if (!hy %in% inputYears) {
    stop("harmonizationPeriod ", hy, " is not a year in the input data")
  }
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  # total area of each cell, which must be constant over time and equal in input and target
  targetArea <- dimSums(setYears(xTarget[, hy, ], NULL), dim = 3)
  stopifnot(all(abs(dimSums(xTarget, dim = 3) - targetArea) < 10^-5))

  # apply absolute changes of input data to target data of the harmonization year
  changed <- setYears(xTarget[, hy, ], NULL) +
    (xInput[, inputYears > hy, ] - setYears(xInput[, hy, ], NULL))
  stopifnot(all(abs(dimSums(changed, dim = 3) - targetArea) < 10^-5))

  # absolute changes can become negative if the input data loses more area of a
  # category than the target data has in the harmonization year
  forest <- intersect(c("pltns", "primf", "secdf"), getItems(changed, dim = 3))
  otherLand <- intersect(c("primn", "secdn"), getItems(changed, dim = 3))

  if (length(forest) > 0 && any(changed[, , forest] < 0)) {
    negativeValues <- as.vector(changed[, , forest])
    negativeValues <- negativeValues[negativeValues < 0]
    toolStatusMessage("warn", paste0("absolute changes made forest categories negative: ",
                                     length(negativeValues), " of ", length(changed), " cells (",
                                     round(100 * length(negativeValues) / length(changed), 1),
                                     "%), min = ", min(negativeValues),
                                     ", mean = ", mean(negativeValues),
                                     ", median = ", median(negativeValues),
                                     ", shortfall was deducted from ",
                                     paste(otherLand, collapse = "/")))
    # deduct missing area from other land categories instead of letting forest
    # categories become negative, proportional to their current shares
    if (length(otherLand) > 0) {
      shortfall <- -changed[, , forest] * (changed[, , forest] < 0)
      available <- dimSums(changed[, , otherLand], dim = 3)
      totalShortfall <- dimSums(shortfall, dim = 3)
      deducted <- pmin(totalShortfall, available)
      changed[, , forest] <- changed[, , forest] +
        shortfall * deducted / (totalShortfall + (totalShortfall == 0))
      changed[, , otherLand] <- changed[, , otherLand] -
        deducted * changed[, , otherLand] / (available + (available == 0))
    }
  }

  if (any(changed < 0)) {
    changed[changed < 0] <- 0
    prim <- intersect(c("primf", "primn"), getItems(changed, dim = 3))
    urban <- intersect("urban", getItems(changed, dim = 3))
    other <- setdiff(getItems(changed, dim = 3), c(prim, urban))
    sumOf <- function(items) {
      if (length(items) > 0) {
        dimSums(changed[, , items], dim = 3)
      } else {
        0
      }
    }
    primSum <- sumOf(prim)
    urbanSum <- sumOf(urban)
    otherSum <- sumOf(other)
    factor <- (targetArea - primSum - urbanSum) / (otherSum + (otherSum == 0))
    factor[!is.finite(factor) | factor < 0] <- 0
    if (length(other) > 0) {
      changed[, , other] <- changed[, , other] * pmin(factor, 1)
    }
    if (any(primSum + urbanSum > targetArea + 10^-5)) {
      toolStatusMessage("warn", paste0("prim + urban area exceed total area after correcting ",
                                       "negative values, scaling prim categories down"))
    }
    factor <- (targetArea - urbanSum) / (primSum + (primSum == 0))
    factor[!is.finite(factor) | factor < 0] <- 0
    if (length(prim) > 0) {
      changed[, , prim] <- changed[, , prim] * pmin(factor, 1)
    }
    if (length(urban) > 0 && any(urbanSum > targetArea + 10^-5)) {
      stop("urban area exceeds total area after correcting negative values, ",
           "urban would need to be scaled but that is not allowed")
    }
  }

  out <- mbind(xTarget[, targetYears <= hy, ], changed)

  # during harmonization primf and primn expansion might be introduced due to
  # primf or primn differences between input and target dataset
  # replace primf and primn expansion with secdf and secdn
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

  stopifnot(all(abs(dimSums(out, dim = 3) - targetArea) < 10^-5))

  return(out)
}
