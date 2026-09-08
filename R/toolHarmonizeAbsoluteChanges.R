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
#' If this makes values negative, they are set to 0 and the other categories
#' are reduced proportionally, so that the total area stays constant.
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
  out <- mbind(xTarget[, targetYears <= hp, ],
               setYears(xTarget[, hp, ], NULL) +
                 (xInput[, inputYears > hp, ] - setYears(xInput[, hp, ], NULL)))

  # absolute changes can become negative if the input data loses more area of a
  # category than the target data has in the harmonization year, set these to 0
  # and reduce the remaining categories proportionally so that the total area
  # remains unchanged
  if (any(out < 0)) {
    toolStatusMessage("warn", paste0("negative values introduced by absolute changes were set to 0, ",
                                     "remaining categories were reduced proportionally to keep the ",
                                     "total area constant"),
                      level = 1)
    out[out < 0] <- 0
    correction <- dimSums(setYears(xTarget[, hp, ], NULL), dim = 3) / dimSums(out, dim = 3)
    correction[!is.finite(correction)] <- 1
    out <- out * correction
  }

  # during harmonization primf and primn expansion might be introduced due to
  # primf or primn differences between input and target dataset
  # replace primf and primn expansion with secdf and secdn
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100, level = 3)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100, level = 3)

  return(out)
}
