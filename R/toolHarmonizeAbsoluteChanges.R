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
#' set to 0 and the non-prim categories are then scaled, so that the total area
#' stays constant.
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
  hp <- harmonizationPeriod

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)

  stopifnot(length(hp) == 1,
            round(hp) == hp,
            setequal(getItems(xInput, dim = 1), getItems(xTarget, dim = 1)),
            setequal(getItems(xInput, dim = 3), getItems(xTarget, dim = 3)))
  if (!hp %in% targetYears) {
    stop("harmonizationPeriod ", hp, " is not a year in the target data")
  }
  if (!hp %in% inputYears) {
    stop("harmonizationPeriod ", hp, " is not a year in the input data")
  }
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  # apply absolute changes of input data to target data of the harmonization year
  changed <- setYears(xTarget[, hp, ], NULL) +
    (xInput[, inputYears > hp, ] - setYears(xInput[, hp, ], NULL))

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
                                     paste(otherLand, collapse = "/")),
                      level = 1)
    # deduct missing area from other land categories instead of letting forest
    # categories become negative, proportional to their current shares
    for (from in forest) {
      shortfall <- -changed[, , from] * (changed[, , from] < 0)
      if (length(otherLand) > 0) {
        available <- dimSums(changed[, , otherLand], dim = 3)
        deducted <- pmin(shortfall, available)
        changed[, , from] <- changed[, , from] + deducted
        safeAvailable <- available + (available == 0)
        for (to in otherLand) {
          changed[, , to] <- changed[, , to] - deducted * changed[, , to] / safeAvailable
        }
      }
    }
  }

  # set remaining negative values to 0 and scale non-prim categories afterwards
  # so that the total area remains unchanged (primf and primn are left untouched)
  if (any(changed < 0)) {
    changed[changed < 0] <- 0
    targetArea <- dimSums(setYears(xTarget[, hp, ], NULL), dim = 3)
    prim <- intersect(c("primf", "primn"), getItems(changed, dim = 3))
    nonPrim <- setdiff(getItems(changed, dim = 3), c("primf", "primn"))
    primSum <- if (length(prim) > 0) dimSums(changed[, , prim], dim = 3) else 0
    factor <- (targetArea - primSum) / dimSums(changed[, , nonPrim], dim = 3)
    if (any(factor < 0, na.rm = TRUE)) {
      toolStatusMessage("warn", paste0("prim area exceeds total area after correcting negative ",
                                       "values, setting non-prim categories to 0 in affected cells"),
                        level = 1)
    }
    factor[!is.finite(factor) | factor < 0] <- 0
    changed[, , nonPrim] <- changed[, , nonPrim] * factor
    changed[is.na(changed)] <- 0
  }

  out <- mbind(xTarget[, targetYears <= hp, ], changed)

  # during harmonization primf and primn expansion might be introduced due to
  # primf or primn differences between input and target dataset
  # replace primf and primn expansion with secdf and secdn
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100, level = 3)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100, level = 3)

  return(out)
}
