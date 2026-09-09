#' toolHarmonizeAbsoluteChanges
#'
#' Harmonize datasets by applying the absolute changes of the input data relative
#' to a single harmonization year on top of the target value at that year:
#' out[t] = target[harmonizationYear] + (input[t] - input[harmonizationYear]) for
#' t >= harmonizationYear and out[t] = target[t] for t < harmonizationYear.
#' This is the aneris/mip "offset" method without convergence (a permanent offset).
#' If forest categories would become negative, the shortfall is funded from
#' primn/secdn, remaining negatives are clamped to zero and then non-prim items
#' are scaled so that the total area is conserved.
#'
#' @param xInput input data as magpie object
#' @param xTarget target data as magpie object
#' @param harmonizationPeriod single integer year, the harmonization year,
#' which has to exist in both input and target data
#' @param level passed to toolStatusMessage and toolReplaceExpansion
#' @return harmonized data set as magpie object with target data for years
#' before the harmonization year and the target value at the harmonization year
#' plus the absolute changes of the input data for later years
#' @author Pascal Sauer
toolHarmonizeAbsoluteChanges <- function(xInput, xTarget, harmonizationPeriod, level = 3) {
  y <- harmonizationPeriod

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)
  futureYears <- inputYears[inputYears > y]

  stopifnot(length(y) == 1,
            round(y) == y,
            y %in% targetYears,
            y %in% inputYears,
            futureYears %in% inputYears,
            setequal(getItems(xInput, dim = 1), getItems(xTarget, dim = 1)),
            setequal(getItems(xInput, dim = 3), getItems(xTarget, dim = 3)))
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  # apply absolute changes of input data relative to the harmonization year
  # on top of the target value at the harmonization year
  delta <- xInput[, futureYears, ] - setYears(xInput[, y, ], NULL)
  changed <- setYears(xTarget[, y, ], NULL) + delta
  out <- mbind(xTarget[, targetYears <= y, ], changed)

  # negative handling
  forest <- intersect(c("pltns", "primf", "secdf"), getItems(out, dim = 3))
  otherLand <- intersect(c("primn", "secdn"), getItems(out, dim = 3))
  negForest <- as.vector(out[, , forest])
  negForest <- negForest[negForest < 0]

  if (length(negForest) > 0) {
    toolStatusMessage("warn", paste0(length(negForest), " negative forest cells (",
                                     round(100 * length(negForest) / length(as.vector(out)), 1),
                                     " % of all cells), negative values: min = ", signif(min(negForest), 3),
                                     ", mean = ", signif(mean(negForest), 3),
                                     ", median = ", signif(median(negForest), 3)),
                      level = level)
    # fund the shortfall from primn/secdn instead of letting forest go negative,
    # deducting proportional to their current shares
    future <- out[, futureYears, ]
    for (f in forest) {
      shortfall <- future[, , f]
      shortfall[shortfall > 0] <- 0
      shortfall <- -shortfall
      available <- dimSums(future[, , otherLand], dim = 3)
      deduction <- pmin(shortfall, available)
      # avoid 0/0 = NaN for cells where deduction and available are both 0
      available[available == 0] <- 1
      future[, , otherLand] <- future[, , otherLand] - deduction * future[, , otherLand] / available
      future[, , f] <- future[, , f] + deduction
    }
    out[, futureYears, ] <- future
  }

  # clamp remaining negatives, this raises affected totals
  out[out < 0] <- 0

  # restore totals by scaling non-prim items, leaving primf/primn untouched
  prim <- intersect(c("primf", "primn"), getItems(out, dim = 3))
  nonPrim <- setdiff(getItems(out, dim = 3), c("primf", "primn"))
  targetArea <- setYears(dimSums(xTarget[, y, ], dim = 3), NULL)
  future <- out[, futureYears, ]
  primSum <- if (length(prim) > 0) dimSums(future[, , prim], dim = 3) else 0
  nonPrimSum <- dimSums(future[, , nonPrim], dim = 3)
  factor <- (targetArea - primSum) / nonPrimSum
  if (any(as.vector(factor) < 0, na.rm = TRUE)) {
    toolStatusMessage("warn", "prim areas exceed total area, non-prim areas are set to 0",
                      level = level)
    factor[factor < 0] <- 0
  }
  future[, , nonPrim] <- future[, , nonPrim] * factor
  future[is.na(future)] <- 0
  out[, futureYears, ] <- future

  # negative handling can introduce primf/primn expansion, replace it here
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100, level = level)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100, level = level)

  return(out)
}
