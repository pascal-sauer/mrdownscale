test_that("toolRescaleRegionTotals restores the expected region totals", {
  cells <- c("regone.1", "regone.2", "regtwo.1")
  xAmount <- new.magpie(cells, years = 2025, names = "amount.c3ann", fill = 0)
  xAmount["regone.1", , ] <- 60
  xAmount["regtwo.1", , ] <- 10
  expectedRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  expectedRegionAmount["regone", , ] <- 70
  expectedRegionAmount["regtwo", , ] <- 10
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(expectedRegionAmount) <- c("region", "year", "category", "data")

  conditions <- captureConditions(rescaled <- toolRescaleRegionTotals(xAmount, expectedRegionAmount))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # regone lost amounts on its zero cell and rescales its own remaining cell by
  # 7/6, regtwo already holds its expected total and stays untouched
  expect_equal(as.vector(rescaled["regone.1", , ]), 70)
  expect_equal(as.vector(rescaled["regone.2", , ]), 0)
  expect_equal(as.vector(rescaled["regtwo.1", , ]), 10)
  expect_equal(as.numeric(dimSums(rescaled, 1.2)), as.numeric(expectedRegionAmount))
})

test_that("toolRescaleRegionTotals falls back to a global factor scoped to the offending table", {
  cells <- c("regone.1", "regtwo.1")
  xAmount <- new.magpie(cells, years = c(2025, 2030), names = "amount.c3ann", fill = 0)
  xAmount["regone.1", 2030, ] <- 10
  xAmount["regtwo.1", , ] <- 10
  expectedRegionAmount <- new.magpie(c("regone", "regtwo"), years = c(2025, 2030),
                                     names = "amount.c3ann", fill = 0)
  expectedRegionAmount["regone", c(2025, 2030), ] <- c(10, 15)
  expectedRegionAmount["regtwo", , ] <- 10
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(expectedRegionAmount) <- c("region", "year", "category", "data")

  conditions <- captureConditions(rescaled <- toolRescaleRegionTotals(xAmount, expectedRegionAmount))
  # regone expects 10 units in 2025 but holds none, so the warn check fires
  expect_true(any(grepl("Check failed.*hosting capacity", conditions$conditions)))

  # the 2025 table falls back to the global factor 20/10, doubling regtwo while
  # regone keeps holding nothing
  expect_equal(as.vector(rescaled["regone.1", 2025, ]), 0)
  expect_equal(as.vector(rescaled["regtwo.1", 2025, ]), 20)

  # the 2030 table has no unholdable region and stays on its regional factor 15/10
  expect_equal(as.vector(rescaled["regone.1", 2030, ]), 15)
  expect_equal(as.vector(rescaled["regtwo.1", 2030, ]), 10)
  expect_equal(as.vector(dimSums(rescaled, 1.2)["regone", 2030, ]), 15)
})
