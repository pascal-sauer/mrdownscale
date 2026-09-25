items <- c("primf", "primn", "secdf", "secdn", "urban", "other")

test_that("toolHarmonizeAbsoluteChanges works", {
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

  # after the harmonization year absolute changes from input are applied:
  # secdf grows by 2 Mha from 2020 to 2025 -> target secdf 10 Mha + 2 = 12 Mha
  expect_equal(as.vector(out["reg.one", 2025, ]), c(40, 10, 12, 5, 5, 28))
  expect_equal(as.vector(out["reg.one", 2030, ]), c(39, 10, 14, 5, 5, 27))

  # primf expansion in the input data is replaced with secdf
  expect_equal(as.vector(out["reg.two", 2025, ]), c(30, 10, 23, 5, 5, 27))
  expect_equal(as.vector(out["reg.two", 2030, ]), c(30, 10, 25, 5, 5, 25))

  # total area is constant over time and unchanged
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 8))

  # invalid harmonizationPeriod
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2015),
               "hy %in% inputYears is not TRUE")
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2010, 2020)))
})

test_that("toolHarmonizeAbsoluteChanges avoids negative values", {
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
  xTarget <- new.magpie("reg.four", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.four", year, ] <- c(40, 10, 10, 5, 5, 30)
  }

  xInput <- new.magpie("reg.four", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.four", 2020, ] <- c(20, 10, 20, 5, 10, 35)
  # input loses 15 Mha secdf, target only has 10 Mha secdf in 2020
  xInput["reg.four", 2025, ] <- c(20, 10, 5, 5, 10, 50)

  out <- expect_message(
    toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020),
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

test_that("toolHarmonizeAbsoluteChanges scales remaining categories while protecting prim and urban", {
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
  # and only the remaining categories are scaled down to keep the total area constant
  expect_equal(as.vector(out["reg.five", 2025, c("primf", "primn", "secdf", "secdn")]), c(40, 0, 0, 0))
  expect_equal(as.vector(out["reg.five", 2025, c("urban", "other")]), c(5, 55))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges scales prim if prim and urban exceed the total area", {
  xTarget <- new.magpie("reg.seven", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.seven", year, ] <- c(40, 4, 10, 1, 5, 40)
  }

  xInput <- new.magpie("reg.seven", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.seven", 2020, ] <- c(20, 4, 20, 1, 5, 50)
  # prim gains of 76 Mha make prim overshoot the total area once negative
  # categories are clipped to 0
  xInput["reg.seven", 2025, ] <- c(90, 10, 0, 0, 0, 0)

  # the forest notification is a message and the prim notification a warning, so
  # run the harmonizer once per expected condition
  out <- expect_message(
    suppressWarnings(
      toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
    ),
    "absolute changes made forest categories negative"
  )
  expect_warning(
    suppressMessages(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)),
    "prim \\+ urban exceed total area"
  )

  # prim is scaled down from 110 to 100, toolReplaceExpansion afterwards moves the
  # primf expansion of 60 Mha into secdf
  expect_equal(as.vector(out["reg.seven", 2025, ]), c(40, 0, 60, 0, 0, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges refuses to scale urban", {
  xTarget <- new.magpie("reg.eight", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.eight", year, ] <- c(0, 0, 10, 5, 80, 5)
  }

  xInput <- new.magpie("reg.eight", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.eight", 2020, ] <- c(0, 0, 10, 5, 10, 75)
  # urban gain of 90 Mha makes urban alone overshoot the total area once other
  # land is clipped to 0, which cannot be fixed without scaling urban
  xInput["reg.eight", 2025, ] <- c(0, 0, 0, 0, 100, 0)

  expect_error(
    suppressWarnings(suppressMessages(
      toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
    )),
    "urban would need to be scaled but that is not allowed"
  )
})

test_that("toolHarmonizeAbsoluteChanges rejects inconsistent or invalid input data", {
  xTarget <- new.magpie("reg.six", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.six", year, ] <- c(40, 10, 10, 5, 5, 30)
  }

  xInput <- new.magpie("reg.six", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.six", 2020, ] <- c(20, 10, 20, 5, 10, 35)
  xInput["reg.six", 2025, ] <- c(20, 10, 22, 5, 10, 33)

  # total area of input is not constant over time
  brokenInput <- xInput
  brokenInput["reg.six", 2025, "other"] <- 30
  expect_error(toolHarmonizeAbsoluteChanges(brokenInput, xTarget, 2020))

  # total area of target is not constant over time
  brokenTarget <- xTarget
  brokenTarget["reg.six", 2010, "other"] <- 25
  expect_error(toolHarmonizeAbsoluteChanges(xInput, brokenTarget, 2020))

  # input data contains NA values
  naInput <- xInput
  naInput["reg.six", 2025, "urban"] <- NA_real_
  expect_error(toolHarmonizeAbsoluteChanges(naInput, xTarget, 2020))
})

test_that("toolGetHarmonizer returns the absoluteChanges harmonizer", {
  harmonizer <- toolGetHarmonizer("absoluteChanges")
  expect_true(is.function(harmonizer))
  expect_error(toolGetHarmonizer("nonexistent"))
})
