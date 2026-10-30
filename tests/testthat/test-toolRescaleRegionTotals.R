test_that("toolRescaleRegionTotals restores the target region totals", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  xAmount <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  xAmount["regone.1", , ] <- 60
  xAmount["regtwo.1", , ] <- 10
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 70
  targetRegionAmount["regtwo", , ] <- 10
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")

  conditions <- captureConditions(rescaled <- toolRescaleRegionTotals(xAmount, targetRegionAmount))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # regone lost amounts on its zero cluster and rescales its own remaining cluster
  # by 7/6, regtwo already holds its target total and stays untouched
  expect_equal(as.vector(rescaled["regone.1", , ]), 70)
  expect_equal(as.vector(rescaled["regone.2", , ]), 0)
  expect_equal(as.vector(rescaled["regtwo.1", , ]), 10)
  expect_equal(as.numeric(dimSums(rescaled, 1.2)), as.numeric(targetRegionAmount))
})

test_that("toolRescaleRegionTotals falls back to a global factor scoped to the offending variable", {
  clusters <- c("regone.1", "regtwo.1")
  xAmount <- new.magpie(clusters, years = c(2025, 2030), names = "amount.c3ann", fill = 0)
  xAmount["regone.1", 2030, ] <- 10
  xAmount["regtwo.1", , ] <- 10
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = c(2025, 2030),
                                   names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", c(2025, 2030), ] <- c(10, 15)
  targetRegionAmount["regtwo", , ] <- 10
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")

  conditions <- captureConditions(rescaled <- toolRescaleRegionTotals(xAmount, targetRegionAmount))
  # regone expects 10 units in 2025 but holds none, so the warn check fires
  expect_true(any(grepl("Check failed.*receiving capacity", conditions$conditions)))

  # 2025 falls back to the global factor 20/10, doubling regtwo while
  # regone keeps holding nothing
  expect_equal(as.vector(rescaled["regone.1", 2025, ]), 0)
  expect_equal(as.vector(rescaled["regtwo.1", 2025, ]), 20)

  # 2030 has no unholdable region and stays on its regional factor 15/10
  expect_equal(as.vector(rescaled["regone.1", 2030, ]), 15)
  expect_equal(as.vector(rescaled["regtwo.1", 2030, ]), 10)
  expect_equal(as.vector(dimSums(rescaled, 1.2)["regone", 2030, ]), 15)
})
