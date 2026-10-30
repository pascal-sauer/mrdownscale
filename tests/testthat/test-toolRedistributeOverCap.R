test_that("toolRedistributeOverCap moves the excess into the headroom of the region", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  xAmount <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  xAmount["regone.1", , ] <- 150
  xAmount["regone.2", , ] <- 10
  xAmount["regtwo.1", , ] <- 1
  capAmount <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 120)
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(capAmount) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(redistributed <- toolRedistributeOverCap(xAmount, capAmount))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # regone's 30 units of excess fit into the 110 units of headroom of its own
  # region, so regtwo stays untouched
  expect_equal(as.vector(redistributed["regone.1", , ]), 120)
  expect_equal(as.vector(redistributed["regone.2", , ]), 40)
  expect_equal(as.vector(redistributed["regtwo.1", , ]), 1)
  expect_equal(as.numeric(dimSums(redistributed, c(1, 3))), as.numeric(dimSums(xAmount, c(1, 3))))
})

test_that("toolRedistributeOverCap falls back to the global pool with a warning", {
  clusters <- c("regone.1", "regtwo.1")
  xAmount <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  xAmount["regone.1", , ] <- 150
  xAmount["regtwo.1", , ] <- 10
  capAmount <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 120)
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(capAmount) <- c("region", "id", "year", "category", "data")

  # regone has no headroom of its own, its 30 units of excess fit into regtwo's
  # 110 units of spare over the global pool, reported as a warning
  conditions <- captureConditions(redistributed <- toolRedistributeOverCap(xAmount, capAmount))
  expect_true(any(grepl("Check failed.*Regions absorb their own excess", conditions$conditions)))
  expect_equal(as.vector(redistributed["regone.1", , ]), 120)
  expect_equal(as.vector(redistributed["regtwo.1", , ]), 40)
})

test_that("toolRedistributeOverCap leaves infeasible variables unchanged", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  xAmount <- new.magpie(clusters, years = c(2025, 2030), names = "amount.c3ann", fill = 0)
  xAmount["regone.1", , ] <- 150
  xAmount["regone.2", , ] <- 10
  xAmount["regtwo.1", , ] <- 10
  capAmount <- new.magpie(clusters, years = c(2025, 2030), names = "amount.c3ann", fill = 120)
  capAmount["regone.2", 2025, ] <- 10
  capAmount["regtwo.1", , ] <- 12
  getSets(xAmount) <- c("region", "id", "year", "category", "data")
  getSets(capAmount) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(redistributed <- toolRedistributeOverCap(xAmount, capAmount))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # in 2025 not even the global pool has enough headroom for the 30 units of
  # excess, so it stays unchanged instead of losing amounts
  expect_equal(as.vector(redistributed["regone.1", 2025, ]), 150)
  expect_equal(as.vector(redistributed["regone.2", 2025, ]), 10)
  expect_equal(as.vector(redistributed["regtwo.1", 2025, ]), 10)

  # in 2030 the excess fits into the headroom of its own region
  expect_equal(as.vector(redistributed["regone.1", 2030, ]), 120)
  expect_equal(as.vector(redistributed["regone.2", 2030, ]), 40)
  expect_equal(as.vector(redistributed["regtwo.1", 2030, ]), 10)
})
