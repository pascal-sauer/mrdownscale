#' calcLandHarmonized
#'
#' This function computes a version of the chosen land input data in the resolution
#' of the land input data but harmonized to the land use information of the
#' chosen land target data set (harmonized categories as well as harmonized
#' transition from historic target data to simulated input data).
#'
#' @inheritParams calcLandInput
#' @param target name of the land target source to be used
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets.
#' For harmonization = "absoluteChanges" this must instead be a single
#' integer value, the harmonization year, which must be present in the target
#' dataset.
#' @param harmonization name of harmonization method, see \code{\link{toolGetHarmonizer}}
#' @author Pascal Sauer, Jan Philipp Dietrich
calcLandHarmonized <- function(input, target, harmonizationPeriod, harmonization) {
  xInput <- calcOutput("LandInputRecategorized", input = input, target = target, aggregate = FALSE)
  geometry <- attr(xInput, "geometry")
  crs <- attr(xInput, "crs")

  if ((harmonization == "absoluteChanges") != (length(harmonizationPeriod) == 1)) {
    stop("harmonizationPeriod must be a single year for harmonization = \"absoluteChanges\" ",
         "and a vector of two years for the other harmonization methods")
  }

  # absoluteChanges only uses the target data up to the harmonization year,
  # no extrapolation is needed
  if (harmonization == "absoluteChanges") {
    xTarget <- calcOutput("LandTargetLowRes", input = input, target = target,
                          endOfHistory = harmonizationPeriod, aggregate = FALSE)
  } else {
    xTarget <- calcOutput("LandTargetExtrapolated", input = input, target = target,
                          harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
  }

  # checks and corrections
  inSum <- dimSums(xInput, dim = 3)
  tSum <- dimSums(xTarget, dim = 3)
  toolExpectLessDiff(inSum, inSum[, 1, ], 10^-5, "Total areas in input stay constant over time")
  toolExpectLessDiff(tSum, tSum[, 1, ], 10^-5, "Total areas in target stay constant over time")
  toolExpectLessDiff(inSum[, 1, ], tSum[, 1, ], 10^-5,
                     "Total areas are the same in target and input data")

  xInput <- toolEqualizeArea(xInput, xTarget[, harmonizationPeriod[1], ])

  harmonizer <- toolGetHarmonizer(harmonization)
  out <- harmonizer(xInput, xTarget, harmonizationPeriod = harmonizationPeriod)

  attr(out, "geometry") <- geometry
  attr(out, "crs")      <- crs

  # checks
  toolExpectTrue(!is.null(attr(out, "geometry")), "Data contains geometry information")
  toolExpectTrue(!is.null(attr(out, "crs")), "Data contains CRS information")
  toolExpectTrue(identical(unname(getSets(out)), c("region", "id", "year", "data")),
                 "Dimensions are named correctly")
  toolExpectTrue(setequal(getItems(out, dim = 3), getItems(xInput, dim = 3)),
                 "Land categories remain unchanged")
  toolExpectTrue(all(out >= 0), "All values are >= 0")
  outSum <- dimSums(out, dim = 3)
  toolExpectLessDiff(outSum, outSum[, 1, ], 10^-5, "Total areas in output stay constant over time")
  toolExpectLessDiff(outSum, dimSums(xInput, dim = 3), 10^-5, "Total areas remain unchanged")
  toolPrimExpansionCheck(out)
  toolExpectLessDiff(out[, getYears(out, as.integer = TRUE) <= harmonizationPeriod[1], ],
                     xTarget[, getYears(xTarget, as.integer = TRUE) <= harmonizationPeriod[1], ],
                     10^-5, "Returning reference data before harmonization period")

  if (length(harmonizationPeriod) == 2) {
    outAfterHarmonization <- out[, getYears(out, as.integer = TRUE) >= harmonizationPeriod[2], ]
    inputAfterHarmonization <- xInput[, getYears(xInput, as.integer = TRUE) >= harmonizationPeriod[2], ]
    nonprimfix <- setdiff(getItems(out, dim = 3), c("primf", "primn", "secdf", "secdn"))
    toolExpectLessDiff(outAfterHarmonization[, , nonprimfix],
                       inputAfterHarmonization[, , nonprimfix],
                       10^-5,
                       "Returning input data after harmonization period (not checking primf/primn/secdf/secdn)")
    toolExpectLessDiff(dimSums(outAfterHarmonization[, , c("primf", "secdf")], 3),
                       dimSums(inputAfterHarmonization[, , c("primf", "secdf")], 3),
                       10^-5, "Returning input data after harmonization period (checking primf + secdf)")
    if ("primn" %in% getItems(out, 3)) {
      toolExpectLessDiff(dimSums(outAfterHarmonization[, , c("primn", "secdn")], 3),
                         dimSums(inputAfterHarmonization[, , c("primn", "secdn")], 3),
                         10^-5,
                         "Returning input data after harmonization period (checking primn + secdn)")
    }
  }

  return(list(x = out,
              class = "magpie",
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Harmonized land data",
              clean_magpie = FALSE))
}
