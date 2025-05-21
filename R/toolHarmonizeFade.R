#' toolHarmonizeFade
#'
#' Tool function for creating a harmonized data set with a smooth s-shaped
#' transition from historic target data to simulated input data.
#'
#' @param xInput input data as magpie object
#' @param xTarget target data as magpie object
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets
#' @return harmonized data set as magpie object with data from input for years
#' before the harmonization period, data from target for years after the
#' harmonization period and a smooth transition in between.
#' @author Jan Philipp Dietrich, Pascal Sauer
toolHarmonizeFade <- function(xInput, xTarget, harmonizationPeriod) {
  a <- harmonizationPeriod[1]
  b <- harmonizationPeriod[2]

  xInput <- toolEqualizeArea(xInput, xTarget)

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)
  transitionYears <- inputYears[inputYears > a & inputYears < b]
  stopifnot(length(harmonizationPeriod) == 2,
            round(harmonizationPeriod) == harmonizationPeriod,
            b > a,
            a %in% targetYears,
            b %in% inputYears,
            transitionYears %in% inputYears,
            transitionYears %in% targetYears,
            getItems(xInput, dim = 1) == getItems(xTarget, dim = 1),
            getItems(xInput, dim = 3) == getItems(xTarget, dim = 3))

  # fade over from extrapolated target data to input data
  out <- convergence(xTarget[, transitionYears, ], xInput[, transitionYears, ],
                     start_year = a, end_year = b, type = "s")
  outYears <- getYears(out, as.integer = TRUE)

  out <- mbind(xTarget[, targetYears < min(outYears), ],
               out,
               xInput[, inputYears > max(outYears), ])

  # during harmonization primf and primn expansion might be introduced due to
  # primf or primn differences between input and target dataset
  # replace primf and primn expansion with secdf and secdn
  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100)

  return(out)
}
