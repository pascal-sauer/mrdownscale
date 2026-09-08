#' toolHarmonizeAbsoluteChanges
#'
#' Harmonize datasets by applying the absolute changes of the input data
#' relative to a single harmonization year on top of the target data at that
#' year. This is the aneris/mip "offset" method without convergence, i.e. a
#' permanent offset that is never faded out.
#'
#' @param xInput magpie object with the data to take absolute changes from,
#' usually model projections
#' @param xTarget magpie object with the historical target data
#' @param harmonizationPeriod a single integer year that must exist in both
#' xInput and xTarget; before this year the target dataset is used, from this
#' year onward the target value at harmonizationPeriod plus the absolute
#' change of the input data (relative to harmonizationPeriod) is used
#' @param level passed to toolReplaceExpansion
#' @return magpie object with harmonized data, negative values are clamped to
#' 0 and each region/timestep is renormalized to keep the total area constant,
#' emitting a warning if any clamping happened
#' @author Pascal Sauer
toolHarmonizeAbsoluteChanges <- function(xInput, xTarget, harmonizationPeriod, level = 1) {
  y <- harmonizationPeriod

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)

  stopifnot(length(y) == 1,
            round(y) == y,
            y %in% targetYears,
            y %in% inputYears,
            setequal(getItems(xInput, dim = 1), getItems(xTarget, dim = 1)),
            setequal(getItems(xInput, dim = 3), getItems(xTarget, dim = 3)))
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  futureYears <- inputYears[inputYears > y]
  delta <- xInput[, futureYears, ] - setYears(xInput[, y, ], NULL)
  changed <- setYears(xTarget[, y, ], NULL) + delta

  out <- mbind(xTarget[, targetYears <= y, ], changed)

  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100, level = level)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100, level = level)

  if (any(out < 0)) {
    toolStatusMessage("warn", paste0("clamping negative values to 0 (min value: ", signif(min(out), 3), ")"),
                      level = level)
    out[out < 0] <- 0
    targetArea <- setYears(dimSums(xTarget[, y, ], dim = 3), NULL)
    out <- out * targetArea / dimSums(out, dim = 3)
    out[is.na(out)] <- 0
  }

  return(out)
}
