items <- c("primf", "primn", "secdf", "secdn", "urban", "pltns",
           "pastr", "range", "c3ann_rainfed", "c4ann_rainfed")

test_that("calcLandHarmonized works with harmonization = absoluteChanges", {
  harmonizationYear <- 2020

  # real MAgPIE/LUH3 source data is not available in CI, so the source data
  # functions are mocked with synthetic data
  xTarget <- new.magpie(c("reg.one", "reg.two", "reg.three"), years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.one", year, ] <- c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0)
    xTarget["reg.two", year, ] <- c(30, 10, 20, 5, 5, 0, 30, 0, 0, 0)
    xTarget["reg.three", year, ] <- c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0)
  }

  xInput <- new.magpie(c("reg.one", "reg.two", "reg.three"),
                       years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput["reg.one", 2020, ] <- c(20, 10, 20, 5, 10, 0, 35, 0, 0, 0)
  xInput["reg.one", 2025, ] <- c(20, 10, 22, 5, 10, 0, 33, 0, 0, 0)
  xInput["reg.one", 2030, ] <- c(19, 10, 24, 5, 10, 0, 32, 0, 0, 0)
  xInput["reg.two", 2020, ] <- c(30, 8, 10, 7, 5, 0, 40, 0, 0, 0)
  xInput["reg.two", 2025, ] <- c(32, 8, 11, 7, 5, 0, 37, 0, 0, 0)
  xInput["reg.two", 2030, ] <- c(33, 8, 12, 7, 5, 0, 35, 0, 0, 0)
  # reg.three has a non-increasing prim input trend but a negative raw secdf in
  # 2025, so the compensation cut on primf/primn makes their raw 2030 recovery
  # an expansion even though the input trend never expanded
  xInput["reg.three", 2020, ] <- c(30, 10, 20, 5, 5, 0, 30, 0, 0, 0)
  xInput["reg.three", 2025, ] <- c(28, 10, 2, 0, 5, 0, 55, 0, 0, 0)
  xInput["reg.three", 2030, ] <- c(26, 10, 20, 5, 5, 0, 34, 0, 0, 0)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"

  # temporarily replaces calcOutput in the package namespace, so the internal
  # source data calls return the synthetic data above and anything else errors
  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandInputRecategorized = xInput,
             LandTargetLowRes = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  x <- expect_warning(
    calcLandHarmonized(input = "magpie", target = "luh3",
                       harmonizationPeriod = c(harmonizationYear, harmonizationYear),
                       harmonization = "absoluteChanges")$x,
    regexp = "primf is expanding considerably"
  )

  expect_false(is.null(attr(x, "geometry")))
  expect_false(is.null(attr(x, "crs")))

  targetYears <- c(2010, 2020)
  inputYears <- c(2020, 2025, 2030)
  expect_equal(getYears(x, as.integer = TRUE),
               sort(c(targetYears, inputYears[inputYears > harmonizationYear])))

  # before and in the harmonization year target data is returned
  expect_true(max(abs(x[, targetYears, ] - xTarget[, targetYears, ])) < 10^-5)
  expect_true(min(x) >= 0)

  # total area is constant over time and equal to the total area of the target data
  targetArea <- dimSums(setYears(xTarget[, harmonizationYear, ], NULL), dim = 3)
  expect_true(max(abs(dimSums(x, dim = 3) - targetArea)) < 10^-5)

  # after the harmonization year absolute changes were applied, except for
  # primf/primn/secdf/secdn (prim expansion is replaced with secd expansion)
  xInputEqualized <- toolEqualizeArea(xInput, xTarget[, harmonizationYear, ])
  yearsAfter <- inputYears[inputYears > harmonizationYear]
  raw <- setYears(xTarget[, harmonizationYear, ], NULL) +
    (xInputEqualized[, yearsAfter, ] - setYears(xInputEqualized[, harmonizationYear, ], NULL))
  difference <- x[, yearsAfter, ] - raw
  # exclude cells where the correction was active, in these all categories are
  # rescaled, not only the ones which became negative
  difference <- difference * (dimSums(raw < 0, dim = 3) == 0)
  difference <- difference[, , c("primf", "primn", "secdf", "secdn"), invert = TRUE]
  expect_equal(max(abs(difference)), 0)

  # reg.three: the compensation cut scales primf/primn down in 2025, and their
  # raw recovery in 2030 is capped and replaced with secdf/secdn expansion
  expect_true(max(abs(x["reg.three", 2025, c("primf", "primn")] - c(95 / 3, 25 / 3))) < 10^-5)
  expect_true(max(abs(x["reg.three", 2030, c("primf", "primn", "secdf", "secdn")] -
                        c(95 / 3, 25 / 3, 10 + 36 - 95 / 3, 5 + 10 - 25 / 3))) < 10^-5)
  expect_true(toolMaxExpansion(x["reg.three", , "primf"]) <= 10^-5)
  expect_true(toolMaxExpansion(x["reg.three", , "primn"]) <= 10^-5)
})

test_that("calcLandHarmonized errors if absoluteChanges gets two different years", {
  expect_error(
    calcLandHarmonized(input = "magpie", target = "luh3",
                       harmonizationPeriod = c(2020, 2050),
                       harmonization = "absoluteChanges"),
    regexp = "harmonizationPeriod must be two equal integer values"
  )
})

test_that("calcLandHarmonized errors if harmonizationPeriod is not two integers", {
  expect_error(
    calcLandHarmonized(input = "magpie", target = "luh3",
                       harmonizationPeriod = 2020,
                       harmonization = "absoluteChanges"),
    regexp = "harmonizationPeriod must always be two integer values"
  )
  expect_error(
    calcLandHarmonized(input = "magpie", target = "luh3",
                       harmonizationPeriod = 2020,
                       harmonization = "fade"),
    regexp = "harmonizationPeriod must always be two integer values"
  )
})
