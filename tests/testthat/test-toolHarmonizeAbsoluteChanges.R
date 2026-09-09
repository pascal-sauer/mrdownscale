# fixture shared by the "core formula" tests: constant crop/past categories,
# no clamping needed
makeBaseFixture <- function(regions = "GLO",
                            targetYears = c(2000, 2010, 2020),
                            targetValues = list(crop = c(40, 45, 50), past = c(60, 55, 50)),
                            inputYears = c(2020, 2025, 2030),
                            inputValues = list(crop = c(55, 58, 53), past = c(45, 42, 47))) {
  categories <- names(targetValues)
  stopifnot(identical(categories, names(inputValues)))

  xTarget <- new.magpie(regions, targetYears, categories, sets = c("region", "year", "data"))
  for (category in categories) {
    xTarget[, targetYears, category] <- targetValues[[category]]
  }

  xInput <- new.magpie(regions, inputYears, categories, sets = c("region", "year", "data"))
  for (category in categories) {
    xInput[, inputYears, category] <- inputValues[[category]]
  }

  list(xTarget = xTarget, xInput = xInput)
}

# fixture with 5 categories where the future timestep needs clamping,
# used to test the redistribution scheme (fix 2)
makeClampFixture <- function(regions = "GLO") {
  categories <- c("primf", "primn", "secdf", "secdn", "crop")
  xTarget <- new.magpie(regions, c(2000, 2020), categories, sets = c("region", "year", "data"))
  xTarget[, 2000, ] <- c(20, 5, 30, 10, 35)
  xTarget[, 2020, ] <- c(20, 5, 30, 10, 35)

  xInput <- new.magpie(regions, c(2020, 2030), categories, sets = c("region", "year", "data"))
  xInput[, 2020, ] <- c(20, 5, 30, 10, 35)
  xInput[, 2030, ] <- c(20, 5, 50, 35, -10) # crop would go negative, needs clamping

  list(xTarget = xTarget, xInput = xInput)
}

# fixture reproducing the fix 1 bug: clamping in one timestep followed by an
# unclamped timestep
makeReorderingFixture <- function(regions = "GLO") {
  categories <- c("primf", "secdf", "crop")
  xTarget <- new.magpie(regions, c(2000, 2020), categories, sets = c("region", "year", "data"))
  xTarget[, 2000, ] <- c(30, 60, 10)
  xTarget[, 2020, ] <- c(30, 60, 10)

  xInput <- new.magpie(regions, c(2020, 2030, 2040), categories, sets = c("region", "year", "data"))
  xInput[, 2020, ] <- c(30, 40, 30)
  xInput[, 2030, ] <- c(30, 70, 0) # -> crop goes negative -> clamping
  xInput[, 2040, ] <- c(30, 40, 30) # -> back to a value that needs no clamping

  list(xTarget = xTarget, xInput = xInput)
}

test_that("toolHarmonizeAbsoluteChanges applies the core formula", {
  f <- makeBaseFixture()
  xTarget <- f$xTarget
  xInput <- f$xInput

  out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)

  expect_equal(getYears(out, as.integer = TRUE), c(2000, 2010, 2020, 2025, 2030))
  expect_identical(as.vector(out[, 2020, ]), as.vector(xTarget[, 2020, ]))
  expect_identical(as.vector(out[, 2025, "crop"]),
                   as.vector(xTarget[, 2020, "crop"] + (xInput[, 2025, "crop"] - xInput[, 2020, "crop"])))
  expect_identical(as.vector(out[, 2025, "past"]),
                   as.vector(xTarget[, 2020, "past"] + (xInput[, 2025, "past"] - xInput[, 2020, "past"])))
  expect_identical(as.vector(out[, 2030, "crop"]),
                   as.vector(xTarget[, 2020, "crop"] + (xInput[, 2030, "crop"] - xInput[, 2020, "crop"])))
})

test_that("toolHarmonizeAbsoluteChanges conserves the total area, including when clamping is triggered", {
  f <- makeClampFixture()

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = 2020)
  }, "clamping")

  targetTotal <- as.vector(dimSums(f$xTarget[, 2020, ], dim = 3))
  for (year in getYears(out, as.integer = TRUE)) {
    expect_equal(as.vector(dimSums(out[, year, ], dim = 3)), targetTotal)
  }
})

test_that("toolHarmonizeAbsoluteChanges passes reference data through before the harmonization year", {
  f <- makeBaseFixture()
  out <- toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = 2020)
  expect_identical(as.vector(out[, c(2000, 2010), ]), as.vector(f$xTarget[, c(2000, 2010), ]))
})

test_that("toolHarmonizeAbsoluteChanges clamps negative values, warns with number and share, keeps total", {
  xTarget <- new.magpie("GLO", c(2000, 2020), c("crop", "past"),
                        sets = c("region", "year", "data"))
  xTarget[, 2000, ] <- c(40, 60)
  xTarget[, 2020, ] <- c(5, 95)

  xInput <- new.magpie("GLO", c(2020, 2030), c("crop", "past"),
                       sets = c("region", "year", "data"))
  xInput[, 2020, ] <- c(50, 50)
  xInput[, 2030, ] <- c(0, 100)

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }, "clamping 1 negative values to 0 \\(50% of the values after the harmonization year, min value: -45\\)")

  expect_true(all(out >= 0))
  targetTotal <- as.vector(dimSums(xTarget[, 2020, ], dim = 3))
  expect_equal(as.vector(dimSums(out[, 2030, ], dim = 3)), targetTotal)
})

test_that("toolHarmonizeAbsoluteChanges primf does not re-expand after a clamped timestep (regression, fix 1)", {
  f <- makeReorderingFixture()

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = 2020)
  }, "clamping")

  expect_equal(as.vector(out[, 2020, "primf"]), 30)
  expect_equal(as.vector(out[, 2030, "primf"]), 30)
  expect_equal(as.vector(out[, 2040, "primf"]), 30)
  expect_lte(toolMaxExpansion(out[, , "primf"]), 0)
})

test_that("toolHarmonizeAbsoluteChanges redistributes clamping surplus only among eligible categories", {
  f <- makeClampFixture()

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = 2020)
  }, "clamping")

  # primf/primn keep exactly their absolute-change values, they are not scaled
  expect_equal(as.vector(out[, 2030, "primf"]), 20)
  expect_equal(as.vector(out[, 2030, "primn"]), 5)

  # surplus is taken from secdf/secdn proportionally to their (post-clamping) values
  preClampSecdf <- 50
  preClampSecdn <- 35
  expect_equal(as.vector(out[, 2030, "secdf"]) / as.vector(out[, 2030, "secdn"]),
               preClampSecdf / preClampSecdn)

  # crop was clamped to 0 and has no headroom, stays 0
  expect_equal(as.vector(out[, 2030, "crop"]), 0)

  # total is conserved
  expect_equal(as.vector(dimSums(out[, 2030, ], dim = 3)),
               as.vector(dimSums(f$xTarget[, 2020, ], dim = 3)))

  # all values are >= 0
  expect_true(all(out >= 0))
})

test_that("toolHarmonizeAbsoluteChanges leaves other regions bit-for-bit unaffected when only one region clamps", {
  clampFixture <- makeBaseFixture(regions = "clampMe",
                                  targetYears = c(2000, 2020),
                                  targetValues = list(crop = c(40, 5), past = c(60, 95)),
                                  inputYears = c(2020, 2030),
                                  inputValues = list(crop = c(50, 0), past = c(50, 100)))
  fineFixture <- makeBaseFixture(regions = "fineRegion",
                                 targetYears = c(2000, 2020),
                                 targetValues = list(crop = c(40, 50), past = c(60, 50)),
                                 inputYears = c(2020, 2030),
                                 inputValues = list(crop = c(55, 53), past = c(45, 47)))

  xTarget <- mbind(clampFixture$xTarget, fineFixture$xTarget)
  xInput <- mbind(clampFixture$xInput, fineFixture$xInput)

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)
  }, "clamping")

  # the region that needed no clamping is returned bit-for-bit unchanged
  expectedFine <- fineFixture$xTarget[, 2020, ] +
    (fineFixture$xInput[, 2030, ] - fineFixture$xInput[, 2020, ])
  expect_identical(as.vector(out["fineRegion", 2030, ]), as.vector(expectedFine))

  # each region's total is conserved separately
  expect_equal(as.vector(dimSums(out["clampMe", 2030, ], dim = 3)),
               as.vector(dimSums(xTarget["clampMe", 2020, ], dim = 3)))
  expect_equal(as.vector(dimSums(out["fineRegion", 2030, ], dim = 3)),
               as.vector(dimSums(xTarget["fineRegion", 2020, ], dim = 3)))

  expect_true(all(out >= 0))
})

test_that("toolHarmonizeAbsoluteChanges returns an unclamped timestep unchanged even if another timestep is clamped", {
  f <- makeReorderingFixture()

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = 2020)
  }, "clamping")

  # 2040 needed no clamping, so it must be exactly the raw absolute-change value (factor exactly 1)
  expected2040 <- f$xTarget[, 2020, ] + (f$xInput[, 2040, ] - f$xInput[, 2020, ])
  expect_identical(as.vector(out[, 2040, ]), as.vector(expected2040))
})

test_that("toolHarmonizeAbsoluteChanges validates its inputs", {
  f <- makeBaseFixture()

  expect_error(toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = 2021),
               "targetYears")

  expect_error(toolHarmonizeAbsoluteChanges(f$xInput, f$xTarget, harmonizationPeriod = c(2020, 2030)),
               "length\\(y\\) == 1")

  noFutureInput <- f$xInput[, 2020, ]
  expect_error(toolHarmonizeAbsoluteChanges(noFutureInput, f$xTarget, harmonizationPeriod = 2020),
               "futureYears")

  inputWithNA <- f$xInput
  inputWithNA[1, 1, 1] <- NA
  expect_error(toolHarmonizeAbsoluteChanges(inputWithNA, f$xTarget, harmonizationPeriod = 2020),
               "anyNA")
})
