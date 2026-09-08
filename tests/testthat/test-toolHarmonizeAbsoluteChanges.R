test_that("toolHarmonizeAbsoluteChanges applies the core formula", {
  xTarget <- new.magpie("GLO", c(2000, 2010, 2020), c("crop", "past"),
                        sets = c("region", "year", "data"))
  xTarget[, 2000, ] <- c(40, 60)
  xTarget[, 2010, ] <- c(45, 55)
  xTarget[, 2020, ] <- c(50, 50)

  xInput <- new.magpie("GLO", c(2020, 2025, 2030), c("crop", "past"),
                       sets = c("region", "year", "data"))
  xInput[, 2020, ] <- c(55, 45)
  xInput[, 2025, ] <- c(58, 42)
  xInput[, 2030, ] <- c(53, 47)

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

test_that("toolHarmonizeAbsoluteChanges conserves the total area", {
  xTarget <- new.magpie("GLO", c(2000, 2010, 2020), c("crop", "past"),
                        sets = c("region", "year", "data"))
  xTarget[, 2000, ] <- c(40, 60)
  xTarget[, 2010, ] <- c(45, 55)
  xTarget[, 2020, ] <- c(50, 50)

  xInput <- new.magpie("GLO", c(2020, 2025, 2030), c("crop", "past"),
                       sets = c("region", "year", "data"))
  xInput[, 2020, ] <- c(55, 45)
  xInput[, 2025, ] <- c(58, 42)
  xInput[, 2030, ] <- c(53, 47)

  out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)

  targetTotal <- as.vector(dimSums(xTarget[, 2020, ], dim = 3))
  for (year in getYears(out, as.integer = TRUE)) {
    expect_equal(as.vector(dimSums(out[, year, ], dim = 3)), targetTotal)
  }
})

test_that("toolHarmonizeAbsoluteChanges passes reference data through before the harmonization year", {
  xTarget <- new.magpie("GLO", c(2000, 2010, 2020), c("crop", "past"),
                        sets = c("region", "year", "data"))
  xTarget[, 2000, ] <- c(40, 60)
  xTarget[, 2010, ] <- c(45, 55)
  xTarget[, 2020, ] <- c(50, 50)

  xInput <- new.magpie("GLO", c(2020, 2025, 2030), c("crop", "past"),
                       sets = c("region", "year", "data"))
  xInput[, 2020, ] <- c(55, 45)
  xInput[, 2025, ] <- c(58, 42)
  xInput[, 2030, ] <- c(53, 47)

  out <- toolHarmonizeAbsoluteChanges(xInput, xTarget, harmonizationPeriod = 2020)

  expect_identical(as.vector(out[, c(2000, 2010), ]), as.vector(xTarget[, c(2000, 2010), ]))
})

test_that("toolHarmonizeAbsoluteChanges clamps negative values and warns, keeping the total conserved", {
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
  }, "clamping negative values")

  expect_true(all(out >= 0))
  targetTotal <- as.vector(dimSums(xTarget[, 2020, ], dim = 3))
  expect_equal(as.vector(dimSums(out[, 2030, ], dim = 3)), targetTotal)
})
