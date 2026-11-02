test_that("toolAllocateRegionTotals distributes targets over clusters proportional to their capacity", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 30
  targetRegionAmount["regtwo", , ] <- 12
  clusterCapacity <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterCapacity["regone.1", , ] <- 60
  clusterCapacity["regone.2", , ] <- 20
  clusterCapacity["regtwo.1", , ] <- 24
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterCapacity) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterCapacity))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # every region reaches exactly its target, its clusters share it proportional
  # to their capacity (regone 30 over capacities 60:20, regtwo 12 on its one cluster)
  expect_equal(as.vector(allocated["regone.1", , ]), 22.5)
  expect_equal(as.vector(allocated["regone.2", , ]), 7.5)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 12)
  expect_equal(as.vector(dimSums(allocated, 1.2)), as.vector(targetRegionAmount))
})

test_that("toolAllocateRegionTotals pools the target of regions without capacity globally", {
  clusters <- c("regone.1", "regone.2", "regtwo.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 20
  targetRegionAmount["regtwo", , ] <- 10
  clusterCapacity <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterCapacity["regone.1", , ] <- 60
  clusterCapacity["regone.2", , ] <- 30
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterCapacity) <- c("region", "id", "year", "category", "data")

  # regtwo targets 10 but has no capacity at all to take it, its target is
  # pooled globally over the 90 capacity of regone (60:30 on top of 20 own)
  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterCapacity))
  expect_true(any(grepl("Check failed.*Regions with target totals > 0 have receiving capacity",
                        conditions$conditions)))
  expect_equal(as.vector(allocated["regone.1", , ]), (20 + 10) * 60 / 90)
  expect_equal(as.vector(allocated["regone.2", , ]), (20 + 10) * 30 / 90)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 0)
  expect_equal(as.numeric(dimSums(allocated, c(1, 3))), as.numeric(dimSums(targetRegionAmount, c(1, 3))))
})

test_that("toolAllocateRegionTotals spreads stranded targets only over the remaining room", {
  clusters <- c("regone.1", "regone.2", "regtwo.1", "regthree.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo", "regthree"), years = 2025,
                                   names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 90
  targetRegionAmount["regtwo", , ] <- 10
  targetRegionAmount["regthree", , ] <- 20
  clusterCapacity <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterCapacity["regone.1", , ] <- 60
  clusterCapacity["regone.2", , ] <- 30
  clusterCapacity["regthree.1", , ] <- 40
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterCapacity) <- c("region", "id", "year", "category", "data")

  # regone's own target fills both of its clusters exactly, so the stranded 10
  # of regtwo without capacity may only go to the remaining room at regthree
  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterCapacity))
  expect_true(any(grepl("Check failed.*Regions with target totals > 0 have receiving capacity",
                        conditions$conditions)))
  expect_equal(as.vector(allocated["regone.1", , ]), 60)
  expect_equal(as.vector(allocated["regone.2", , ]), 30)
  expect_equal(as.vector(allocated["regtwo.1", , ]), 0)
  expect_equal(as.vector(allocated["regthree.1", , ]), 30)
  expect_true(all(as.vector(allocated) <= as.vector(clusterCapacity) + 10^-10))
})

test_that("toolAllocateRegionTotals has nowhere to put stranded targets where all capacities are zero", {
  clusters <- c("regone.1", "regtwo.1")
  targetRegionAmount <- new.magpie(c("regone", "regtwo"), years = 2025, names = "amount.c3ann", fill = 0)
  targetRegionAmount["regone", , ] <- 5
  targetRegionAmount["regtwo", , ] <- 7
  clusterCapacity <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterCapacity) <- c("region", "id", "year", "category", "data")

  # without any global room the pooled target cannot be allocated at all (the
  # caller has to account for the loss), and must not turn into NaN
  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterCapacity))
  expect_true(any(grepl("Check failed.*Regions with target totals > 0 have receiving capacity",
                        conditions$conditions)))
  expect_equal(as.vector(allocated), c(0, 0))
  expect_false(anyNA(as.vector(allocated)))
})

test_that("toolAllocateRegionTotals never puts anything on clusters without capacity", {
  clusters <- c("regone.1", "regone.2")
  targetRegionAmount <- new.magpie("regone", years = 2025, names = "amount.c3ann", fill = 24)
  clusterCapacity <- new.magpie(clusters, years = 2025, names = "amount.c3ann", fill = 0)
  clusterCapacity["regone.1", , ] <- 8
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterCapacity) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterCapacity))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # the whole 24 of the region goes to the only cluster with capacity
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
  clusterCapacity <- new.magpie(clusters, years = c(2025, 2030),
                                names = c("amount.c3ann", "amount.c4ann"), fill = 0)
  clusterCapacity["regone.1", , ] <- 9
  clusterCapacity["regone.2", , ] <- 3
  clusterCapacity["regtwo.1", , "amount.c3ann"] <- 16
  clusterCapacity["regone", 2030, "amount.c4ann"] <- 0
  getSets(targetRegionAmount) <- c("region", "year", "category", "data")
  getSets(clusterCapacity) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(allocated <- toolAllocateRegionTotals(targetRegionAmount, clusterCapacity))
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # each variable is allocated on its own: within a variable the clusters share
  # the regional target proportional to their capacities of exactly that variable
  expect_equal(as.vector(allocated["regone.1", 2025, "amount.c3ann"]), 9)
  expect_equal(as.vector(allocated["regone.2", 2025, "amount.c3ann"]), 3)
  expect_equal(as.vector(allocated["regone.1", 2030, "amount.c3ann"]), 22.5)
  expect_equal(as.vector(allocated["regone.1", 2025, "amount.c4ann"]), 4.5)
  expect_equal(as.vector(allocated["regone.2", 2025, "amount.c4ann"]), 1.5)
  # regone targets 6 of c4ann 2030 but the region has no c4ann capacity there,
  # its target is stranded in a pool without any capacity and vanishes
  expect_equal(as.vector(allocated["regone.1", 2030, "amount.c4ann"]), 0)
  expect_equal(as.vector(allocated["regtwo.1", , "amount.c3ann"]), c(8, 8))
  expect_equal(as.vector(allocated["regtwo.1", , "amount.c4ann"]), c(0, 0))
  expect_true(min(allocated) >= 0)
})
