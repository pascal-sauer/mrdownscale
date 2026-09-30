items <- c("primf", "primn", "secdf", "secdn", "urban", "other",
           "pltns", "pastr", "range", "c3ann_rainfed", "c4ann_rainfed")

test_that("toolHarmonizeAbsoluteChanges works", {
  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.one", year, ] <- c(40, 10, 10, 5, 5, 30, 0, 0, 0, 0, 0)
    xTarget["reg.two", year, ] <- c(30, 10, 20, 5, 5, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput["reg.one", 2020, ] <- c(20, 10, 20, 5, 10, 35, 0, 0, 0, 0, 0)
  xInput["reg.one", 2025, ] <- c(20, 10, 22, 5, 10, 33, 0, 0, 0, 0, 0)
  xInput["reg.one", 2030, ] <- c(19, 10, 24, 5, 10, 32, 0, 0, 0, 0, 0)
  xInput["reg.two", 2020, ] <- c(30, 8, 10, 7, 5, 40, 0, 0, 0, 0, 0)
  xInput["reg.two", 2025, ] <- c(32, 8, 11, 7, 5, 37, 0, 0, 0, 0, 0)
  xInput["reg.two", 2030, ] <- c(33, 8, 12, 7, 5, 35, 0, 0, 0, 0, 0)

  suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  })

  expect_equal(getYears(out, as.integer = TRUE), c(2010, 2020, 2025, 2030))

  # before and at the harmonization year target data is used
  expect_equal(as.vector(out["reg.one", 2010, ]), c(40, 10, 10, 5, 5, 30, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.two", 2010, ]), c(30, 10, 20, 5, 5, 30, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.one", 2020, ]), c(40, 10, 10, 5, 5, 30, 0, 0, 0, 0, 0))

  # after the harmonization year absolute changes from input are applied:
  # secdf grows by 2 Mha from 2020 to 2025 -> target secdf 10 Mha + 2 = 12 Mha
  expect_equal(as.vector(out["reg.one", 2025, ]), c(40, 10, 12, 5, 5, 28, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.one", 2030, ]), c(39, 10, 14, 5, 5, 27, 0, 0, 0, 0, 0))

  # primf expansion in the input data is replaced with secdf
  expect_equal(as.vector(out["reg.two", 2025, ]), c(30, 10, 23, 5, 5, 27, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.two", 2030, ]), c(30, 10, 25, 5, 5, 25, 0, 0, 0, 0, 0))

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
    xTarget["reg.three", year, ] <- c(5, 5, 10, 5, 5, 70, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.three", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.three", 2020, ] <- c(50, 5, 10, 5, 5, 25, 0, 0, 0, 0, 0)
  # input loses 10 Mha primf, but target only has 5 Mha primf in 2020
  xInput["reg.three", 2025, ] <- c(40, 5, 10, 5, 5, 35, 0, 0, 0, 0, 0)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  expect_equal(as.vector(out["reg.three", 2025, "primf"]), 0)
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
  # target data before the harmonization year is untouched
  expect_equal(as.vector(out["reg.three", 2010, ]), c(5, 5, 10, 5, 5, 70, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.three", 2020, ]), c(5, 5, 10, 5, 5, 70, 0, 0, 0, 0, 0))
})

test_that("toolHarmonizeAbsoluteChanges compensates negative forest area within the forest group", {
  xTarget <- new.magpie("reg.four", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.four", year, ] <- c(40, 10, 10, 5, 5, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.four", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.four", 2020, ] <- c(20, 10, 20, 5, 10, 35, 0, 0, 0, 0, 0)
  # input loses 15 Mha secdf, target only has 10 Mha secdf in 2020
  xInput["reg.four", 2025, ] <- c(20, 10, 5, 5, 10, 50, 0, 0, 0, 0, 0)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  # the secdf shortfall of 5 Mha is compensated by scaling down the other
  # categories of the forest group (primf, primn, secdf, secdn) proportionally,
  # secdf itself is set to 0
  expect_equal(as.vector(out["reg.four", 2025, c("primf", "primn", "secdf", "secdn", "urban", "other")]),
               c(400 / 11, 100 / 11, 0, 50 / 11, 5, 45))
  # the forest group keeps its total area
  expect_equal(as.vector(dimSums(out["reg.four", 2025, c("primf", "primn", "secdf", "secdn")], dim = 3)), 50)
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
  expect_equal(as.vector(out["reg.four", 2020, ]), c(40, 10, 10, 5, 5, 30, 0, 0, 0, 0, 0))
})

test_that("toolHarmonizeAbsoluteChanges scales all categories except urban", {
  xTarget <- new.magpie("reg.five", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.five", year, ] <- c(40, 4, 10, 1, 5, 40, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.five", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.five", 2020, ] <- c(20, 4, 20, 1, 5, 50, 0, 0, 0, 0, 0)
  # input loses 20 Mha secdf, more than secdf + primn + secdn can cover within
  # the forest group after clipping
  xInput["reg.five", 2025, ] <- c(20, 4, 0, 1, 5, 70, 0, 0, 0, 0, 0)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  # the forest group is scaled down so that it keeps its total area of 35 Mha,
  # urban is untouched and other land does not need any further scaling
  expect_equal(as.vector(out["reg.five", 2025, c("primf", "primn", "secdf", "secdn")]),
               c(280 / 9, 28 / 9, 0, 7 / 9))
  expect_equal(as.vector(out["reg.five", 2025, c("urban", "other")]), c(5, 60))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges scales down excess area and replaces prim expansion", {
  xTarget <- new.magpie("reg.seven", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.seven", year, ] <- c(40, 4, 10, 1, 5, 40, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.seven", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.seven", 2020, ] <- c(20, 4, 20, 1, 5, 50, 0, 0, 0, 0, 0)
  # prim gains of 76 Mha make prim overshoot the total area once negative
  # categories are clipped to 0
  xInput["reg.seven", 2025, ] <- c(90, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  # after clipping and scaling, toolReplaceExpansion moves the prim expansions
  # into secdf and secdn
  expect_equal(as.vector(out["reg.seven", 2025, ]), c(40, 4, 155 / 3, 13 / 3, 0, 0, 0, 0, 0, 0, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges refuses to scale urban", {
  xTarget <- new.magpie("reg.eight", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.eight", year, ] <- c(0, 0, 10, 5, 80, 5, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.eight", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.eight", 2020, ] <- c(0, 0, 10, 5, 10, 75, 0, 0, 0, 0, 0)
  # urban gain of 90 Mha makes urban alone overshoot the total area once other
  # land is clipped to 0, which cannot be fixed without scaling urban
  xInput["reg.eight", 2025, ] <- c(0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0)

  expect_error(
    suppressWarnings(suppressMessages(
      toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
    )),
    "all(nonUrbanTarget >= 0) is not TRUE", fixed = TRUE
  )
})

test_that("toolHarmonizeAbsoluteChanges compensates negative cropland within the cropland group", {
  xTarget <- new.magpie("reg.crop", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.crop", year, ] <- c(5, 0, 0, 0, 5, 40, 0, 2, 1, 30, 17)
  }

  xInput <- new.magpie("reg.crop", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.crop", 2020, ] <- c(5, 0, 0, 0, 5, 20, 0, 2, 1, 50, 17)
  # input loses 40 Mha of c3ann_rainfed, more than the 30 Mha in the target
  xInput["reg.crop", 2025, ] <- c(5, 0, 0, 0, 5, 20, 0, 2, 1, 10, 57)

  suppressWarnings(suppressMessages({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }))

  # the c3ann_rainfed shortfall is covered by c4ann_rainfed, the cropland group
  # keeps its total area of 47 Mha
  expect_equal(as.vector(out["reg.crop", 2025, c("c3ann_rainfed", "c4ann_rainfed")]), c(0, 47))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges rejects inconsistent or invalid input data", {
  xTarget <- new.magpie("reg.six", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.six", year, ] <- c(40, 10, 10, 5, 5, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.six", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.six", 2020, ] <- c(20, 10, 20, 5, 10, 35, 0, 0, 0, 0, 0)
  xInput["reg.six", 2025, ] <- c(20, 10, 22, 5, 10, 33, 0, 0, 0, 0, 0)

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
