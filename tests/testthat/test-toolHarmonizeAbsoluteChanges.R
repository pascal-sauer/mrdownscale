test_that("toolHarmonizeAbsoluteChanges works", {
  items <- c("primf", "primn", "secdf", "secdn", "urban", "other")

  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.one", year, ] <- c(40, 10, 10, 5, 5, 30)
    xTarget["reg.two", year, ] <- c(30, 10, 20, 5, 5, 30)
  }

  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput["reg.one", 2020, ] <- c(20, 10, 20, 5, 10, 35)
  xInput["reg.one", 2025, ] <- c(20, 10, 22, 5, 10, 33)
  xInput["reg.one", 2030, ] <- c(19, 10, 24, 5, 10, 32)
  xInput["reg.two", 2020, ] <- c(30, 8, 10, 7, 5, 40)
  xInput["reg.two", 2025, ] <- c(32, 8, 11, 7, 5, 37)
  xInput["reg.two", 2030, ] <- c(33, 8, 12, 7, 5, 35)

  suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  })

  expect_equal(getYears(out, as.integer = TRUE), c(2010, 2020, 2025, 2030))

  # before and at the harmonization year target data is used
  expect_equal(as.vector(out["reg.one", 2010, ]), c(40, 10, 10, 5, 5, 30))
  expect_equal(as.vector(out["reg.two", 2010, ]), c(30, 10, 20, 5, 5, 30))
  expect_equal(as.vector(out["reg.one", 2020, ]), c(40, 10, 10, 5, 5, 30))

  # after the harmonization year absolute changes from input are applied to the
  # target value in the harmonization year (secdf grows by 2 Mha from 2020 to 2025
  # in the input data, target secdf in 2020 is 10 Mha -> 12 Mha in 2025)
  expect_equal(as.vector(out["reg.one", 2025, ]), c(40, 10, 12, 5, 5, 28))
  expect_equal(as.vector(out["reg.one", 2030, ]), c(39, 10, 14, 5, 5, 27))

  # primf expansion in the input data is replaced with secdf
  expect_equal(as.vector(out["reg.two", 2025, ]), c(30, 10, 23, 5, 5, 27))
  expect_equal(as.vector(out["reg.two", 2030, ]), c(30, 10, 25, 5, 5, 25))

  # total area is constant over time and unchanged
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 8))

  # invalid harmonizationPeriod
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2015),
               "2015")
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2010, 2020)))
})

test_that("toolHarmonizeAbsoluteChanges avoids negative values", {
  items <- c("primf", "primn", "secdf", "secdn", "urban", "other")

  xTarget <- new.magpie("reg.three", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.three", year, ] <- c(5, 5, 10, 5, 5, 70)
  }

  xInput <- new.magpie("reg.three", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.three", 2020, ] <- c(50, 5, 10, 5, 5, 25)
  # input loses 10 Mha primf, but target only has 5 Mha primf in 2020
  xInput["reg.three", 2025, ] <- c(40, 5, 10, 5, 5, 35)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  expect_equal(as.vector(out["reg.three", 2025, "primf"]), 0)
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
  # target data before the harmonization year is untouched
  expect_equal(as.vector(out["reg.three", 2010, ]), c(5, 5, 10, 5, 5, 70))
  expect_equal(as.vector(out["reg.three", 2020, ]), c(5, 5, 10, 5, 5, 70))
})

test_that("toolHarmonizeAbsoluteChanges deducts negative forest area from other land", {
  items <- c("primf", "primn", "secdf", "secdn", "urban", "other")

  xTarget <- new.magpie("reg.four", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.four", year, ] <- c(40, 10, 10, 5, 5, 30)
  }

  xInput <- new.magpie("reg.four", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.four", 2020, ] <- c(20, 10, 20, 5, 10, 35)
  # input loses 15 Mha secdf, target only has 10 Mha secdf in 2020
  xInput["reg.four", 2025, ] <- c(20, 10, 5, 5, 10, 50)

  expect_warning(
    {
      out <- suppressMessages(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020))
    },
    "absolute changes made forest categories negative: 1 of 6 cells \\(16.7%\\), min = -5, mean = -5, median = -5"
  )

  # secdf shortfall of 5 Mha is deducted from primn and secdn proportional to
  # their shares, secdf itself is set to 0 and primf stays untouched
  expect_equal(as.vector(out["reg.four", 2025, c("primf", "primn", "secdf", "secdn", "urban", "other")]),
               c(40, 20 / 3, 0, 10 / 3, 5, 45))
  # primn + secdn are reduced by exactly the shortfall
  expect_equal(as.vector(dimSums(out["reg.four", 2025, c("primn", "secdn")], dim = 3)), 10)
  # no clamping or rescaling was needed, other categories are unchanged
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
  expect_equal(as.vector(out["reg.four", 2020, ]), c(40, 10, 10, 5, 5, 30))
})

test_that("toolHarmonizeAbsoluteChanges clamps and scales non-prim categories if other land is insufficient", {
  items <- c("primf", "primn", "secdf", "secdn", "urban", "other")

  xTarget <- new.magpie("reg.five", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.five", year, ] <- c(40, 4, 10, 1, 5, 40)
  }

  xInput <- new.magpie("reg.five", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.five", 2020, ] <- c(20, 4, 20, 1, 5, 50)
  # input loses 20 Mha secdf, primn + secdn (5 Mha) cannot cover the 10 Mha shortfall
  xInput["reg.five", 2025, ] <- c(20, 4, 0, 1, 5, 70)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  # primn and secdn are fully deducted, the remaining secdf deficit is clamped to 0
  # and the non-prim categories are scaled down to keep the total area constant
  expect_equal(as.vector(out["reg.five", 2025, c("primf", "primn", "secdf", "secdn")]), c(40, 0, 0, 0))
  # primf stays exactly at the target value, it is not affected by the scaling
  expect_equal(as.vector(out["reg.five", 2025, c("urban", "other")]), c(60 / 13, 720 / 13))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolGetHarmonizer returns the absoluteChanges harmonizer", {
  harmonizer <- toolGetHarmonizer("absoluteChanges")
  expect_true(is.function(harmonizer))
  expect_error(toolGetHarmonizer("nonexistent"))
})
