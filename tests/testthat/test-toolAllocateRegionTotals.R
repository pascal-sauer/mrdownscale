test_that("toolAllocateRegionTotals spreads targets onto clusters proportional to their weight", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 30
  targetRegionAmount["regtwo", , ] <- 12
  clusterWeight <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterWeight["regone.1", , ] <- 60
  clusterWeight["regone.2", , ] <- 20
  clusterWeight["regtwo.1", , ] <- 24
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterWeight) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterWeight))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # every region holds exactly its target, its clusters share it proportional
  # to their weights (regone 30 over weights 60:20, regtwo 12 on its one cluster)
  expect_equal(as.vector(allocated["regone.1", , ]), 22.5)
  expect_equal(as.vector(allocated["regone.2", , ]), 7.5)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 12)
  expect_equal(as.vector(dimSums(allocated, 1.2)), as.vector(targetRegionAmount))
})

test_that("toolAllocateRegionTotals pools the target of weightless regions globally", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 20
  targetRegionAmount["regtwo", , ] <- 10
  clusterWeight <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterWeight["regone.1", , ] <- 60
  clusterWeight["regone.2", , ] <- 30
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterWeight) <- c("region", "id", "year", "category", "data")

  # regtwo targets 10 but has no weight at all to host it, its target is
  # pooled globally over the 90 weight of regone (60:30 on top of 20 own)
  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterWeight))
  expect_true(any(grepl("Check failed.*Regions with target totals > 0 have hosting capacity",
                        conditions$conditions)))
  expect_equal(as.vector(allocated["regone.1", , ]), (20 + 10) * 60 / 90)
  expect_equal(as.vector(allocated["regone.2", , ]), (20 + 10) * 30 / 90)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 0)
  expect_equal(as.numeric(dimSums(allocated, c(1, 3))), as.numeric(dimSums(targetRegionAmount, c(1, 3))))
})

test_that("toolAllocateRegionTotals never puts anything on weightless clusters", {
  clusters <- c("regone.1", "regone.2")
  targetRegionAmount <- new.magpie("regone", years = 2025, names = "amount.c3ann", fill = 24)
  clusterWeight <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterWeight["regone.1", , ] <- 8
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterWeight) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterWeight))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # the whole 24 of the region falls onto the only weighted cluster
  expect_equal(as.vector(allocated["regone.1", , ]), 24)
  expect_equal(as.vector(allocated["regone.2", , ]), 0)
})

test_that("toolAllocateRegionTotals handles multiple variables and years", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = c(2025, 2030),
                                   names = c("amount.c3ann", "amount.c4ann"), fill = 0)
  targetRegionAmount["regone", 2025, "amount.c3ann"] <- 12
  targetRegionAmount["regone", 2030, "amount.c3ann"] <- 30
  targetRegionAmount["regone", 2025, "amount.c4ann"] <- 6
  targetRegionAmount["regtwo", , "amount.c3ann"] <- 8
  clusterWeight <- new.magpie(clusters, years = c(2025, 2030),
                              names = c("amount.c3ann", "amount.c4ann"), fill = 0)
  clusterWeight["regone.1", , ] <- 9
  clusterWeight["regone.2", , ] <- 3
  clusterWeight["regtwo.1", , "amount.c3ann"] <- 16
  clusterWeight["regone", 2030, "amount.c4ann"] <- 0
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterWeight) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterWeight))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # each variable is allocated on its own: within a variable the clusters share
  # the regional target proportional to their weights of exactly that variable
  expect_equal(as.vector(allocated["regone.1", 2025, "amount.c3ann"]), 9)
  expect_equal(as.vector(allocated["regone.2", 2025, "amount.c3ann"]), 3)
  expect_equal(as.vector(allocated["regone.1", 2030, "amount.c3ann"]), 22.5)
  expect_equal(as.vector(allocated["regone.1", 2025, "amount.c4ann"]), 4.5)
  expect_equal(as.vector(allocated["regone.2", 2025, "amount.c4ann"]), 1.5)
  # regone targets 6 of c4ann 2030 but the region has no c4ann weight there,
  # its target is stranded in a pool without any weight and vanishes
  expect_equal(as.vector(allocated["regone.1", 2030, "amount.c4ann"]), 0)
  expect_equal(as.vector(allocated["regtwo.1", , "amount.c3ann"]), c(8, 8))
  expect_equal(as.vector(allocated["regtwo.1", , "amount.c4ann"]), c(0, 0))
  expect_true(min(allocated) >= 0)
})
