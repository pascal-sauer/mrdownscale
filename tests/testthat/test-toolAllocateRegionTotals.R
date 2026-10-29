test_that("toolAllocateRegionTotals spreads expectations onto cells proportional to their weight", {
  cells <- c("regone.1", "regone.2", "regtwo.1")
  expectedRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  expectedRegionAmount["regone", , ] <- 30
  expectedRegionAmount["regtwo", , ] <- 12
  cellWeight <- new.magpie(cells, years = 2025, names = "amount.c3ann", fill = 0)
  cellWeight["regone.1", , ] <- 60
  cellWeight["regone.2", , ] <- 20
  cellWeight["regtwo.1", , ] <- 24
  getSets(expectedRegionAmount) <- c("region", "year", "category", "data")
  getSets(cellWeight) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(expectedRegionAmount, cellWeight))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # every region holds exactly its expectation, its cells share it proportional
  # to their weights (regone 30 over weights 60:20, regtwo 12 on its one cell)
  expect_equal(as.vector(allocated["regone.1", , ]), 22.5)
  expect_equal(as.vector(allocated["regone.2", , ]), 7.5)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 12)
  expect_equal(as.vector(dimSums(allocated, 1.2)), as.vector(expectedRegionAmount))
})

test_that("toolAllocateRegionTotals pools the expectation of weightless regions globally", {
  cells <- c("regone.1", "regone.2", "regtwo.1")
  expectedRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  expectedRegionAmount["regone", , ] <- 20
  expectedRegionAmount["regtwo", , ] <- 10
  cellWeight <- new.magpie(cells, years = 2025, names = "amount.c3ann", fill = 0)
  cellWeight["regone.1", , ] <- 60
  cellWeight["regone.2", , ] <- 30
  getSets(expectedRegionAmount) <- c("region", "year", "category", "data")
  getSets(cellWeight) <- c("region", "id", "year", "category", "data")

  # regtwo expects 10 but has no weight at all to host it, its expectation is
  # pooled globally over the 90 weight of regone (60:30 on top of 20 own)
  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(expectedRegionAmount, cellWeight))
  expect_true(any(grepl("Check failed.*Regions with expected totals have hosting capacity",
                        conditions$conditions)))
  expect_equal(as.vector(allocated["regone.1", , ]), (20 + 10) * 60 / 90)
  expect_equal(as.vector(allocated["regone.2", , ]), (20 + 10) * 30 / 90)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 0)
  expect_equal(as.numeric(dimSums(allocated, c(1, 3))), as.numeric(dimSums(expectedRegionAmount, c(1, 3))))
})

test_that("toolAllocateRegionTotals never puts anything on weightless cells", {
  cells <- c("regone.1", "regone.2")
  expectedRegionAmount <- new.magpie("regone", years = 2025, names = "amount.c3ann", fill = 24)
  cellWeight <- new.magpie(cells, years = 2025, names = "amount.c3ann", fill = 0)
  cellWeight["regone.1", , ] <- 8
  getSets(expectedRegionAmount) <- c("region", "year", "category", "data")
  getSets(cellWeight) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(expectedRegionAmount, cellWeight))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # the whole 24 of the region falls onto the only weighted cell
  expect_equal(as.vector(allocated["regone.1", , ]), 24)
  expect_equal(as.vector(allocated["regone.2", , ]), 0)
})

test_that("toolAllocateRegionTotals handles multiple items and years per table", {
  cells <- c("regone.1", "regone.2", "regtwo.1")
  expectedRegionAmount <- new.magpie(c("regone", "regtwo"), years = c(2025, 2030),
                                     names = c("amount.c3ann", "amount.c4ann"), fill = 0)
  expectedRegionAmount["regone", 2025, "amount.c3ann"] <- 12
  expectedRegionAmount["regone", 2030, "amount.c3ann"] <- 30
  expectedRegionAmount["regone", 2025, "amount.c4ann"] <- 6
  expectedRegionAmount["regtwo", , "amount.c3ann"] <- 8
  cellWeight <- new.magpie(cells, years = c(2025, 2030),
                           names = c("amount.c3ann", "amount.c4ann"), fill = 0)
  cellWeight["regone.1", , ] <- 9
  cellWeight["regone.2", , ] <- 3
  cellWeight["regtwo.1", , "amount.c3ann"] <- 16
  cellWeight["regone", 2030, "amount.c4ann"] <- 0
  getSets(expectedRegionAmount) <- c("region", "year", "category", "data")
  getSets(cellWeight) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(expectedRegionAmount, cellWeight))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # each table is allocated on its own: within a table the cells share the
  # regional expectation proportional to their weights of exactly that table
  expect_equal(as.vector(allocated["regone.1", 2025, "amount.c3ann"]), 9)
  expect_equal(as.vector(allocated["regone.2", 2025, "amount.c3ann"]), 3)
  expect_equal(as.vector(allocated["regone.1", 2030, "amount.c3ann"]), 22.5)
  expect_equal(as.vector(allocated["regone.1", 2025, "amount.c4ann"]), 4.5)
  expect_equal(as.vector(allocated["regone.2", 2025, "amount.c4ann"]), 1.5)
  # regone expects 6 of c4ann 2030 but the region has no c4ann weight there,
  # its expectation is stranded in a pool without any weight and vanishes
  expect_equal(as.vector(allocated["regone.1", 2030, "amount.c4ann"]), 0)
  expect_equal(as.vector(allocated["regtwo.1", , "amount.c3ann"]), c(8, 8))
  expect_equal(as.vector(allocated["regtwo.1", , "amount.c4ann"]), c(0, 0))
  expect_true(min(allocated) >= 0)
})
