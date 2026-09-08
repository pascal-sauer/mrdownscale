test_that("calcLandHarmonized works with harmonization = absoluteChanges", {
  harmonizationYear <- 2020

  x <- calcOutput("LandHarmonized", input = "magpie", target = "luh3",
                  harmonizationPeriod = harmonizationYear,
                  harmonization = "absoluteChanges", aggregate = FALSE)

  xInput <- calcOutput("LandInputRecategorized", input = "magpie", target = "luh3", aggregate = FALSE)
  xTarget <- calcOutput("LandTargetLowRes", input = "magpie", target = "luh3",
                        endOfHistory = harmonizationYear, aggregate = FALSE)

  expect_false(is.null(attr(x, "geometry")))
  expect_false(is.null(attr(x, "crs")))

  targetYears <- getYears(xTarget, as.integer = TRUE)
  inputYears <- getYears(xInput, as.integer = TRUE)
  expect_equal(getYears(x, as.integer = TRUE),
               sort(c(targetYears[targetYears <= harmonizationYear],
                      inputYears[inputYears > harmonizationYear])))

  # before and in the harmonization year target data is returned
  expect_true(max(x[, targetYears, ] - xTarget) < 10^-5)
  expect_true(min(x) >= 0)

  # total area is constant over time and equal to the total area of the target data
  targetArea <- dimSums(setYears(xTarget[, harmonizationYear, ], NULL), dim = 3)
  expect_true(max(abs(dimSums(x, dim = 3) - targetArea)) < 10^-5)

  # after the harmonization year absolute changes from the input data were applied,
  # except for cells where the raw absolute changes would be negative (these are
  # corrected) and except for primf/primn/secdf/secdn (prim expansion is replaced
  # with secd expansion)
  xInput <- toolEqualizeArea(xInput, xTarget[, harmonizationYear, ])
  yearsAfter <- inputYears[inputYears > harmonizationYear]
  raw <- setYears(xTarget[, harmonizationYear, ], NULL) +
    (xInput[, yearsAfter, ] - setYears(xInput[, harmonizationYear, ], NULL))
  difference <- x[, yearsAfter, ] - raw
  # exclude cells where the correction was active, in these all categories are
  # rescaled, not only the ones which became negative
  difference <- difference * (dimSums(raw < 0, dim = 3) == 0)
  difference <- difference[, , c("primf", "primn", "secdf", "secdn"), invert = TRUE]
  expect_equal(max(abs(difference)), 0)
})

test_that("calcLandHarmonized errors if absoluteChanges gets a harmonization period", {
  expect_error(suppressWarnings(calcOutput(
    "LandHarmonized", input = "magpie", target = "luh3",
    harmonizationPeriod = c(2020, 2050), harmonization = "absoluteChanges", aggregate = FALSE
  )))
})
