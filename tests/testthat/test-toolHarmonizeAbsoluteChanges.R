items <- c("primf", "primn", "secdf", "secdn", "urban", "pltns",
           "pastr", "range", "c3ann_rainfed", "c3ann_irrigated",
           "c4ann_rainfed", "c4ann_irrigated")

# runs the expression and returns its value together with all messages and
# warnings emitted by toolStatusMessage (which uses warnings for "warn" status)
captureConditions <- function(expr) {
  conditions <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      conditions <<- unique(c(conditions, conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      conditions <<- unique(c(conditions, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  return(list(value = value, conditions = conditions))
}

expectCondition <- function(conditions, regexp) {
  expect_true(any(grepl(regexp, conditions)))
}

test_that("toolHarmonizeAbsoluteChanges works", {
  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.one", year, ] <- c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0, 0, 0)
    xTarget["reg.two", year, ] <- c(30, 10, 20, 5, 5, 0, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput["reg.one", 2020, ] <- c(20, 10, 20, 5, 10, 0, 35, 0, 0, 0, 0, 0)
  xInput["reg.one", 2025, ] <- c(20, 10, 22, 5, 10, 0, 33, 0, 0, 0, 0, 0)
  xInput["reg.one", 2030, ] <- c(19, 10, 24, 5, 10, 0, 32, 0, 0, 0, 0, 0)
  xInput["reg.two", 2020, ] <- c(30, 8, 10, 7, 5, 0, 40, 0, 0, 0, 0, 0)
  xInput["reg.two", 2025, ] <- c(32, 8, 11, 7, 5, 0, 37, 0, 0, 0, 0, 0)
  xInput["reg.two", 2030, ] <- c(33, 8, 12, 7, 5, 0, 35, 0, 0, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # no negative absolute changes, so no value is reset to zero and no trend is altered
  # (except prim expansion which is checked below)
  expect_equal(getYears(out, as.integer = TRUE), c(2010, 2020, 2025, 2030))

  # before and at the harmonization year target data is used
  expect_equal(as.vector(out["reg.one", 2010, ]), c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.two", 2010, ]), c(30, 10, 20, 5, 5, 0, 30, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.one", 2020, ]), c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0, 0, 0))

  # after the harmonization year absolute changes from input are applied:
  # secdf grows by 2 Mha from 2020 to 2025 -> target secdf 10 Mha + 2 = 12 Mha
  expect_equal(as.vector(out["reg.one", 2025, ]), c(40, 10, 12, 5, 5, 0, 28, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.one", 2030, ]), c(39, 10, 14, 5, 5, 0, 27, 0, 0, 0, 0, 0))

  # primf expansion in the input data is replaced with secdf
  expect_equal(as.vector(out["reg.two", 2025, ]), c(30, 10, 23, 5, 5, 0, 27, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.two", 2030, ]), c(30, 10, 25, 5, 5, 0, 25, 0, 0, 0, 0, 0))

  # total area is constant over time and unchanged
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 8))

  # invalid harmonizationPeriod
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2015, 2015)),
               "hy %in% inputYears is not TRUE")
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2010, 2020)))
  expect_error(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020))
})

test_that("toolHarmonizeAbsoluteChanges avoids negative values", {
  xTarget <- new.magpie("reg.three", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.three", year, ] <- c(5, 5, 10, 5, 5, 0, 70, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.three", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.three", 2020, ] <- c(50, 5, 10, 5, 5, 0, 25, 0, 0, 0, 0, 0)
  # input loses 10 Mha primf, but target only has 5 Mha primf in 2020
  xInput["reg.three", 2025, ] <- c(40, 5, 10, 5, 5, 0, 35, 0, 0, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  expectCondition(run$conditions, "\\[!\\]")

  expect_equal(as.vector(out["reg.three", 2025, "primf"]), 0)
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
  # target data before the harmonization year is untouched
  expect_equal(as.vector(out["reg.three", 2010, ]), c(5, 5, 10, 5, 5, 0, 70, 0, 0, 0, 0, 0))
  expect_equal(as.vector(out["reg.three", 2020, ]), c(5, 5, 10, 5, 5, 0, 70, 0, 0, 0, 0, 0))
})

test_that("toolHarmonizeAbsoluteChanges compensates negative forest area within the forest group", {
  xTarget <- new.magpie("reg.four", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.four", year, ] <- c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.four", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.four", 2020, ] <- c(20, 10, 20, 5, 10, 0, 35, 0, 0, 0, 0, 0)
  # input loses 15 Mha secdf, target only has 10 Mha secdf in 2020
  xInput["reg.four", 2025, ] <- c(20, 10, 5, 5, 10, 0, 50, 0, 0, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the secdf shortfall of 5 Mha is compensated by scaling down the
  # non-primary categories of the forest group (secdn is scaled to zero),
  # primf and primn stay untouched
  expect_equal(as.vector(out["reg.four", 2025, c("primf", "primn", "secdf", "secdn", "urban", "pastr")]),
               c(40, 10, 0, 0, 5, 45))
  # the forest group keeps its total area
  expect_equal(as.vector(dimSums(out["reg.four", 2025, c("primf", "primn", "secdf", "secdn")], dim = 3)), 50)
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
  expect_equal(as.vector(out["reg.four", 2020, ]), c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0, 0, 0))
})

test_that("toolHarmonizeAbsoluteChanges scales all categories except urban", {
  xTarget <- new.magpie("reg.five", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.five", year, ] <- c(40, 4, 10, 1, 5, 0, 40, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.five", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.five", 2020, ] <- c(20, 4, 20, 1, 5, 0, 50, 0, 0, 0, 0, 0)
  # input loses 20 Mha secdf, more than secdf + primn + secdn can cover within
  # the forest group after clipping
  xInput["reg.five", 2025, ] <- c(20, 4, 0, 1, 5, 0, 70, 0, 0, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the forest group is scaled down so that it keeps its total area of 35 Mha:
  # the non-primary member secdn is scaled to zero first, then primf and primn
  # absorb the remaining 9 Mha proportionally from their 44 Mha total,
  # urban is untouched and pastr does not need any further scaling
  expect_equal(as.vector(out["reg.five", 2025, c("primf", "primn", "secdf", "secdn")]),
               c(40 * 35 / 44, 4 * 35 / 44, 0, 0))
  expect_equal(as.vector(out["reg.five", 2025, c("urban", "pastr")]), c(5, 60))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges protects primary land from compensation cuts", {
  xTarget <- new.magpie("reg.prim", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.prim", year, ] <- c(30, 10, 15, 10, 0, 5, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.prim", years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput["reg.prim", 2020, ] <- c(20, 10, 10, 10, 0, 5, 20, 0, 25, 0, 0, 0)
  # 25 Mha of cropland are converted to pastr, zeroing the orphaned
  # c3ann_rainfed and creating 25 Mha of excess area
  xInput["reg.prim", 2025, ] <- c(20, 10, 10, 10, 0, 5, 45, 0, 0, 0, 0, 0)
  xInput["reg.prim", 2030, ] <- c(20, 10, 10, 10, 0, 5, 20, 0, 25, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the 25 Mha excess is absorbed by secdf, secdn, pltns and pastr (factor 12/17),
  # primf and primn stay at their target values
  expect_equal(as.vector(out["reg.prim", 2025, ]),
               c(30, 10, 180 / 17, 120 / 17, 0, 60 / 17, 660 / 17, 0, 0, 0, 0, 0))
  # since primf was not cut in 2025, its recovery in 2030 is not capped by
  # toolReplaceExpansion and the target data is reached again
  expect_equal(as.vector(out["reg.prim", 2030, ]), c(30, 10, 15, 10, 0, 5, 30, 0, 0, 0, 0, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 4))
  expect_true(!any(grepl("replaced primf expansion", run$conditions)))
})

test_that("toolHarmonizeAbsoluteChanges scales down excess area and replaces prim expansion", {
  xTarget <- new.magpie("reg.seven", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.seven", year, ] <- c(40, 4, 10, 1, 5, 0, 40, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.seven", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.seven", 2020, ] <- c(20, 4, 20, 1, 5, 0, 50, 0, 0, 0, 0, 0)
  # prim gains of 76 Mha make prim overshoot the total area once negative
  # categories are clipped to 0
  xInput["reg.seven", 2025, ] <- c(90, 10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # after clipping and scaling, toolReplaceExpansion moves the prim expansions
  # into secdf and secdn
  expect_equal(as.vector(out["reg.seven", 2025, ]), c(40, 4, 155 / 3, 13 / 3, 0, 0, 0, 0, 0, 0, 0, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges refuses to scale urban", {
  xTarget <- new.magpie("reg.eight", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.eight", year, ] <- c(0, 0, 10, 5, 80, 0, 5, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.eight", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.eight", 2020, ] <- c(0, 0, 10, 5, 10, 0, 75, 0, 0, 0, 0, 0)
  # urban gain of 90 Mha makes urban alone overshoot the total area once pastr
  # is clipped to 0, which cannot be fixed without scaling urban
  xInput["reg.eight", 2025, ] <- c(0, 0, 0, 0, 100, 0, 0, 0, 0, 0, 0, 0)

  expect_error(
    suppressWarnings(suppressMessages(
      toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020))
    )),
    "all(nonUrbanTarget >= -10^-5) is not TRUE", fixed = TRUE
  )
})

test_that("toolHarmonizeAbsoluteChanges sets a group with a negative total area to zero", {
  xTarget <- new.magpie("reg.nine", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.nine", year, ] <- c(50, 0, 0, 0, 0, 0, 50, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.nine", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.nine", 2020, ] <- c(90, 0, 0, 0, 0, 0, 10, 0, 0, 0, 0, 0)
  # primf loses 90 Mha, more than the 50 Mha of the target, so the total area of
  # the forest group becomes negative and the whole group is set to zero
  xInput["reg.nine", 2025, ] <- c(0, 0, 0, 0, 0, 0, 100, 0, 0, 0, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  expect_equal(as.vector(out["reg.nine", 2025, c("primf", "urban", "pastr")]), c(0, 0, 100))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges compensates negative cropland within the cropland group", {
  xTarget <- new.magpie("reg.crop", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.crop", year, ] <- c(5, 0, 0, 0, 5, 0, 42, 1, 30, 0, 17, 0)
  }

  xInput <- new.magpie("reg.crop", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.crop", 2020, ] <- c(5, 0, 0, 0, 5, 0, 22, 1, 50, 0, 17, 0)
  # input loses 40 Mha of c3ann_rainfed, more than the 30 Mha in the target
  xInput["reg.crop", 2025, ] <- c(5, 0, 0, 0, 5, 0, 22, 1, 10, 0, 57, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the c3ann_rainfed shortfall is covered by c4ann_rainfed, the cropland group
  # keeps its total area of 47 Mha
  expect_equal(as.vector(out["reg.crop", 2025, c("c3ann_rainfed", "c4ann_rainfed")]), c(0, 47))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges rejects inconsistent or invalid input data", {
  xTarget <- new.magpie("reg.six", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.six", year, ] <- c(40, 10, 10, 5, 5, 0, 30, 0, 0, 0, 0, 0)
  }

  xInput <- new.magpie("reg.six", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.six", 2020, ] <- c(20, 10, 20, 5, 10, 0, 35, 0, 0, 0, 0, 0)
  xInput["reg.six", 2025, ] <- c(20, 10, 22, 5, 10, 0, 33, 0, 0, 0, 0, 0)

  # total area of input is not constant over time
  brokenInput <- xInput
  brokenInput["reg.six", 2025, "pastr"] <- 30
  expect_error(toolHarmonizeAbsoluteChanges(brokenInput, xTarget, c(2020, 2020)))

  # total area of target is not constant over time
  brokenTarget <- xTarget
  brokenTarget["reg.six", 2010, "pastr"] <- 25
  expect_error(toolHarmonizeAbsoluteChanges(xInput, brokenTarget, c(2020, 2020)))

  # input data contains NA values
  naInput <- xInput
  naInput["reg.six", 2025, "urban"] <- NA_real_
  expect_error(toolHarmonizeAbsoluteChanges(naInput, xTarget, c(2020, 2020)))
})

test_that("toolHarmonizeAbsoluteChanges works with nonland data", {
  nonlandItems <- c("bioh.indu", "bioh.animal")
  # nonland data is not constant-sum over time, this must not be rejected
  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget["reg.one", 2010, ] <- c(1, 5)
  xTarget["reg.one", 2020, ] <- c(10, 2)
  xTarget["reg.two", 2010, ] <- c(3, 4)
  xTarget["reg.two", 2020, ] <- c(8, 7)

  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025, 2030), names = nonlandItems, fill = 0)
  xInput["reg.one", 2020, ] <- c(20, 6)
  # reg.one bioh.animal raw value is 2 + (3 - 6) = -1, so the clamp at zero is exercised
  xInput["reg.one", 2025, ] <- c(25, 3)
  xInput["reg.one", 2030, ] <- c(30, 0)
  xInput["reg.two", 2020, ] <- c(4, 10)
  xInput["reg.two", 2025, ] <- c(6, 12)
  xInput["reg.two", 2030, ] <- c(5, 8)

  expect_no_error(
    out <- suppressWarnings(suppressMessages(
      toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020))
    ))
  )
  expect_equal(getYears(out, as.integer = TRUE), c(2010, 2020, 2025, 2030))

  # until the harmonization year the target is used unchanged
  expect_equal(as.vector(out["reg.one", 2010, ]), c(1, 5))
  expect_equal(as.vector(out["reg.one", 2020, ]), c(10, 2))
  expect_equal(as.vector(out["reg.two", 2020, ]), c(8, 7))

  # afterwards the raw absolute changes of the input are applied, clamped at zero
  expect_equal(as.vector(out["reg.one", 2025, ]), c(15, 0))
  expect_equal(as.vector(out["reg.one", 2030, ]), c(20, 0))
  expect_equal(as.vector(out["reg.two", 2025, ]), c(10, 9))
  expect_equal(as.vector(out["reg.two", 2030, ]), c(9, 5))

  expect_true(all(out >= 0))
})

test_that("toolGetHarmonizer returns the absoluteChanges harmonizer", {
  harmonizer <- toolGetHarmonizer("absoluteChanges")
  expect_true(is.function(harmonizer))
  expect_error(toolGetHarmonizer("nonexistent"))
})

test_that("toolHarmonizeAbsoluteChanges compensates negative crop by its rainfed/irrigated counterpart", {
  xTarget <- new.magpie("reg.comp", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.comp", year, ] <- c(0, 0, 0, 0, 0, 0, 40, 0, 30, 20, 7, 3)
  }

  xInput <- new.magpie("reg.comp", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.comp", 2020, ] <- c(0, 0, 0, 0, 0, 0, 10, 0, 70, 10, 7, 3)
  # input loses 40 Mha of c3ann_rainfed, more than the 30 Mha in the target,
  # the gain goes to c4ann_irrigated
  xInput["reg.comp", 2025, ] <- c(0, 0, 0, 0, 0, 0, 10, 0, 30, 10, 7, 43)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # raw c3ann_rainfed drops to -10 Mha, the deficit is absorbed by
  # c3ann_irrigated, all other categories are untouched
  expect_equal(as.vector(out["reg.comp", 2025, c("pastr", "c3ann_rainfed", "c3ann_irrigated",
                                                 "c4ann_rainfed", "c4ann_irrigated")]),
               c(40, 0, 10, 7, 43))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges scales other crops if counterpart cannot compensate fully", {
  xTarget <- new.magpie("reg.comp2", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.comp2", year, ] <- c(0, 0, 0, 0, 0, 0, 40, 0, 30, 5, 20, 5)
  }

  xInput <- new.magpie("reg.comp2", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.comp2", 2020, ] <- c(0, 0, 0, 0, 0, 0, 15, 0, 55, 5, 20, 5)
  # input loses 50 Mha of c3ann_rainfed, the 5 Mha of c3ann_irrigated are not
  # enough to compensate the 20 Mha raw deficit
  xInput["reg.comp2", 2025, ] <- c(0, 0, 0, 0, 0, 0, 65, 0, 5, 5, 20, 5)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # both c3ann categories are zeroed, the remaining excess is covered by
  # scaling down the other crops, pastr is untouched
  expect_equal(as.vector(out["reg.comp2", 2025, c("pastr", "c3ann_rainfed", "c3ann_irrigated",
                                                  "c4ann_rainfed", "c4ann_irrigated")]),
               c(90, 0, 0, 8, 2))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges zeroes crop pair if both are negative", {
  xTarget <- new.magpie("reg.comp3", years = c(2010, 2020), names = items, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.comp3", year, ] <- c(0, 0, 0, 0, 0, 0, 60, 0, 10, 5, 20, 5)
  }

  xInput <- new.magpie("reg.comp3", years = c(2020, 2025), names = items, fill = 0)
  xInput["reg.comp3", 2020, ] <- c(0, 0, 0, 0, 0, 0, 20, 0, 35, 25, 20, 0)
  # input loses 30 Mha of c3ann_rainfed and 25 Mha of c3ann_irrigated, both raw
  # values become negative, the other crops cover the total deficit of 40 Mha
  xInput["reg.comp3", 2025, ] <- c(0, 0, 0, 0, 0, 0, 20, 0, 5, 0, 45, 30)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  expect_equal(as.vector(out["reg.comp3", 2025, c("pastr", "c3ann_rainfed", "c3ann_irrigated",
                                                  "c4ann_rainfed", "c4ann_irrigated")]),
               c(60, 0, 0, 22.5, 17.5))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

bioItems <- c("primf", "primn", "secdf", "secdn", "urban", "pltns", "pastr", "range",
              "c3ann_rainfed_biofuel_1st_gen", "c3ann_irrigated_biofuel_1st_gen",
              "c4ann_rainfed_biofuel_2nd_gen", "c4ann_irrigated_biofuel_2nd_gen")

test_that("toolHarmonizeAbsoluteChanges compensates negative biofuel crop by its counterpart", {
  xTarget <- new.magpie("reg.bio", years = c(2010, 2020), names = bioItems, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.bio", year, ] <- c(0, 0, 0, 0, 0, 0, 40, 0, 30, 20, 7, 3)
  }

  xInput <- new.magpie("reg.bio", years = c(2020, 2025), names = bioItems, fill = 0)
  xInput["reg.bio", 2020, ] <- c(0, 0, 0, 0, 0, 0, 10, 0, 70, 10, 7, 3)
  # same scenario as the base-crop pair test: c3ann_rainfed_biofuel_1st_gen loses
  # 40 Mha, more than its 20 Mha irrigated twin can cover, c4ann absorbs the rest
  xInput["reg.bio", 2025, ] <- c(0, 0, 0, 0, 0, 0, 10, 0, 30, 10, 7, 43)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the biofuel pair self-compensates exactly like a base-crop pair: c3ann rainfed is
  # zeroed and its irrigated twin shrinks to the pair total, c4ann untouched
  expect_equal(as.vector(out["reg.bio", 2025, c("pastr", "c3ann_rainfed_biofuel_1st_gen",
                                                "c3ann_irrigated_biofuel_1st_gen",
                                                "c4ann_rainfed_biofuel_2nd_gen",
                                                "c4ann_irrigated_biofuel_2nd_gen")]),
               c(40, 0, 10, 7, 43))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges pairs base and biofuel crops of a type separately", {
  colItems <- c("primf", "primn", "secdf", "secdn", "urban", "pltns", "pastr", "range",
                "c3ann_rainfed", "c3ann_irrigated",
                "c3ann_rainfed_biofuel_1st_gen", "c3ann_irrigated_biofuel_1st_gen")
  xTarget <- new.magpie("reg.col", years = c(2010, 2020), names = colItems, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.col", year, ] <- c(0, 0, 0, 0, 0, 0, 71, 0, 8, 10, 6, 5)
  }

  xInput <- new.magpie("reg.col", years = c(2020, 2025), names = colItems, fill = 0)
  xInput["reg.col", 2020, ] <- c(0, 0, 0, 0, 0, 0, 40, 0, 25, 10, 20, 5)
  # both raw base and biofuel c3ann_rainfed go negative, each compensated by its own
  # irrigated twin without leaking into the other class
  xInput["reg.col", 2025, ] <- c(0, 0, 0, 0, 0, 0, 65, 0, 10, 10, 10, 5)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # base pair total 3 Mha, biofuel pair total 1 Mha, each reduced onto its
  # irrigated twin only, pastr untouched
  expect_equal(as.vector(out["reg.col", 2025, c("pastr", "c3ann_rainfed", "c3ann_irrigated",
                                                "c3ann_rainfed_biofuel_1st_gen",
                                                "c3ann_irrigated_biofuel_1st_gen")]),
               c(96, 0, 3, 0, 1))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges harmonizes rainfed/irrigated orphans gracefully", {
  orphanItems <- c("primf", "primn", "secdf", "secdn", "urban", "pltns", "pastr", "range",
                   "c3ann_rainfed")
  xTarget <- new.magpie("reg.orphan", years = c(2010, 2020), names = orphanItems, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.orphan", year, ] <- c(40, 10, 10, 5, 5, 0, 30, 0, 0)
  }

  xInput <- new.magpie("reg.orphan", years = c(2020, 2025), names = orphanItems, fill = 0)
  xInput["reg.orphan", 2020, ] <- c(20, 10, 20, 5, 10, 0, 35, 0, 0)
  xInput["reg.orphan", 2025, ] <- c(20, 10, 22, 5, 10, 0, 33, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the orphaned c3ann_rainfed does not abort harmonization, and since no raw
  # value is negative the absolute changes are reproduced exactly
  expect_equal(as.vector(out["reg.orphan", 2025, ]), c(40, 10, 12, 5, 5, 0, 28, 0, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges compensates a negative orphan within its group", {
  orphanItems <- c("primf", "primn", "secdf", "secdn", "urban", "pltns", "pastr", "range",
                   "c3ann_rainfed")
  xTarget <- new.magpie("reg.orphan2", years = c(2010, 2020), names = orphanItems, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.orphan2", year, ] <- c(40, 10, 10, 5, 5, 0, 25, 0, 5)
  }

  xInput <- new.magpie("reg.orphan2", years = c(2020, 2025), names = orphanItems, fill = 0)
  xInput["reg.orphan2", 2020, ] <- c(20, 10, 20, 5, 10, 0, 25, 0, 10)
  # input loses 10 Mha of the orphaned c3ann_rainfed, more than the 5 Mha in the
  # target, so its raw value becomes -5 Mha
  xInput["reg.orphan2", 2025, ] <- c(20, 10, 20, 5, 10, 0, 35, 0, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the -5 Mha orphan is zeroed within the cropland group and the excess area is
  # taken from all categories except urban, scaling primf and primn only as a last
  # resort (factor 0.9 applied to the non-primary categories, prim untouched)
  expect_equal(as.vector(out["reg.orphan2", 2025, ]), c(40, 10, 9, 4.5, 5, 0, 31.5, 0, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges compensates pairs before passing orphans to the group", {
  mixedItems <- c("primf", "primn", "secdf", "secdn", "urban", "pltns", "pastr", "range",
                  "c3ann_rainfed", "c3ann_irrigated", "c4ann_rainfed")
  xTarget <- new.magpie("reg.orphan3", years = c(2010, 2020), names = mixedItems, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.orphan3", year, ] <- c(0, 0, 0, 0, 0, 0, 47, 0, 30, 20, 3)
  }

  xInput <- new.magpie("reg.orphan3", years = c(2020, 2025), names = mixedItems, fill = 0)
  xInput["reg.orphan3", 2020, ] <- c(0, 0, 0, 0, 0, 0, 10, 0, 70, 10, 10)
  # c3ann_rainfed goes 40 Mha below its input level (raw -10), c4ann_rainfed as
  # an orphan goes 10 Mha below its target (raw -7)
  xInput["reg.orphan3", 2025, ] <- c(0, 0, 0, 0, 0, 0, 60, 0, 30, 10, 0)

  run <- captureConditions(toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020)))
  out <- run$value

  # the complete c3ann pair compensates internally first (total 10 Mha on the
  # irrigated twin), then the orphaned c4ann_rainfed deficit brings the whole
  # cropland group down to its raw total of 3 Mha, pastr untouched
  expect_equal(as.vector(out["reg.orphan3", 2025, c("pastr", "c3ann_rainfed", "c3ann_irrigated",
                                                    "c4ann_rainfed")]),
               c(97, 0, 3, 0))
  expect_true(all(out >= 0))
  expect_equal(as.vector(dimSums(out, dim = 3)), rep(100, 3))
})

test_that("toolHarmonizeAbsoluteChanges rejects cropland names without underscore before water regime", {
  looseItems <- c("primf", "primn", "secdf", "secdn", "urban", "pltns", "pastr", "range",
                  "c3ann_rainfed", "c3ann_irrigated", "rainfed_rest")
  xTarget <- new.magpie("reg.loose", years = c(2010, 2020), names = looseItems, fill = 0)
  for (year in c(2010, 2020)) {
    xTarget["reg.loose", year, ] <- c(40, 10, 10, 5, 5, 0, 20, 0, 5, 5, 0)
  }

  xInput <- new.magpie("reg.loose", years = c(2020, 2025), names = looseItems, fill = 0)
  xInput["reg.loose", 2020, ] <- c(20, 10, 20, 5, 10, 0, 15, 0, 10, 10, 5)
  xInput["reg.loose", 2025, ] <- c(20, 10, 22, 5, 10, 0, 13, 0, 10, 10, 5)

  # rainfed_rest contains the water regime token without the leading underscore,
  # so it belongs to no category group and must be rejected
  expect_error(
    suppressWarnings(suppressMessages(
      toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = c(2020, 2020))
    )),
    "setequal(c(unlist(groups), \"urban\"), getItems(changed, 3)) is not TRUE", fixed = TRUE
  )
})
