harvestItems <- paste0("wood_harvest_area.", c("primf", "secyf", "secmf", "pltns", "primn", "secnf"))
landItems <- c("primf", "secdf", "pltns", "primn", "secdn")

makeHarvest <- function(years, cellValues) {
  x <- new.magpie(names(cellValues), years = years, names = harvestItems, fill = 0)
  getSets(x) <- c("region", "id", "year", "category", "data")
  for (cell in names(cellValues)) {
    for (item in names(cellValues[[cell]])) {
      x[cell, , paste0("wood_harvest_area.", item)] <- cellValues[[cell]][[item]]
    }
  }
  return(x)
}

makeLand <- function(years, cellValues) {
  x <- new.magpie(names(cellValues), years = years, names = landItems, fill = 0)
  getSets(x) <- c("region", "id", "year", "land")
  for (cell in names(cellValues)) {
    for (item in names(cellValues[[cell]])) {
      x[cell, , item] <- cellValues[[cell]][[item]]
    }
  }
  return(x)
}

test_that("absoluteChanges keeps regional and global harvest area parallel to the input", {
  # regone grows its secondary harvest in the input while its harmonized land
  # stays flat, regtwo is inert: the region and global harvest area totals must
  # equal the target total in the harmonization year plus the regional input
  # delta, the prim harvest stays land-derived (1 Mha yr-1 on regone.1) and the
  # secondary harvest absorbs the difference within the region
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = c(100, 100, 95, 90), primn = 40, secdf = 50, secdn = 20, pltns = 10),
    regone.2 = list(primf = 40, primn = 20, secdf = 40, secdn = 10, pltns = 10),
    regtwo.1 = list(primf = 30, primn = 10, secdf = 30, secdn = 15, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 8, 14), pltns = 1, secnf = 1),
                             regone.2 = anchor, regtwo.1 = anchor))

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandHarmonized = land,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )
  conditions <- captureConditions(
    result <- calcWoodHarvestAreaHarmonized(input = "magpie", target = "luh3",
                                            harmonizationPeriod = c(2020, 2020),
                                            harmonization = "absoluteChanges")
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  x <- result$x

  # regional and global totals over all land categories: target total in 2020
  # (8 per region with two cells, 4 in regtwo) plus the input delta
  regionTotal <- dimSums(dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2), 3)
  expect_equal(as.vector(regionTotal["regone", 2025, ]), 8 + (14 - 8))
  expect_equal(as.vector(regionTotal["regone", 2030, ]), 8 + (20 - 8))
  expect_equal(as.vector(regionTotal["regtwo", 2025, ]), 4)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2025, ]), 14 + 4)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2030, ]), 20 + 4)

  # the prim harvest stays exactly land-derived, the 1 Mha yr-1 that the prim
  # harvest already over-achieves scales down both regone secondary cells
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.primf"]), 1, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 8 * 13 / 14, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 2 * 13 / 14, tolerance = 10^-8)
  # inert regtwo keeps its anchor, the target history is untouched
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secyf"]), 2)
  expect_equal(as.vector(x[, 2020, harvestItems]), as.vector(xTarget[, 2020, harvestItems]))
})

test_that("absoluteChanges compensates clamped losses within the same region", {
  # regone.1's raw absolute change turns negative and is clamped to zero, so
  # its region total falls short of the expectation and must be compensated by
  # the other cells of the same region: expected 5.5 over an actual 6 scales
  # regone.2 from 2 to 1.833, while inert regtwo keeps its anchor of 2
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regone.2 = list(primf = 20, primn = 5, secdf = 40, secdn = 10, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 15, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  # regone.1's target is only 1 but its input loses 1.5 from 2 to 0.5, so its
  # raw absolute change of -0.5 turns negative and gets clamped to 0
  xTarget["regone.1", , "wood_harvest_area.secyf"] <- 1
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 0.5, 0.5), pltns = 1, secnf = 1),
                             regone.2 = anchor, regtwo.1 = anchor))

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandHarmonized = land,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )
  conditions <- captureConditions(
    result <- calcWoodHarvestAreaHarmonized(input = "magpie", target = "luh3",
                                            harmonizationPeriod = c(2020, 2020),
                                            harmonization = "absoluteChanges")
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  x <- result$x

  # the clamped cell's secondary harvest stays at zero, its region partner
  # compensates the loss
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 0)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 2 * 5.5 / 6, tolerance = 10^-8)
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secyf"]), 2)

  # the region total stays parallel: 7 + (6.5 - 8) = 5.5
  regionTotal <- dimSums(dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2), 3)
  expect_equal(as.vector(regionTotal["regone", 2025, ]), 5.5, tolerance = 10^-8)
  expect_equal(as.vector(regionTotal["regtwo", 2025, ]), 4)
})

test_that("absoluteChanges falls back to the global pool if a region cannot host the harvest", {
  # regone expects 22 aggregated secondary harvest in 2025 but its secondary
  # land only allows 15, so the excess of regone.1 goes over the global pool
  # into regtwo's headroom: the warn check fires, the region is off-parallel,
  # but the global total stays parallel while the global capacity suffices
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 10, secdn = 5, pltns = 5),
    regone.2 = list(primf = 20, primn = 5, secdf = 5, secdn = 5, pltns = 5),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 5, pltns = 5)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 20, 20), pltns = 1, secnf = 1),
                             regone.2 = anchor, regtwo.1 = anchor))

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandHarmonized = land,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )
  conditions <- captureConditions(
    result <- calcWoodHarvestAreaHarmonized(input = "magpie", target = "luh3",
                                            harmonizationPeriod = c(2020, 2020),
                                            harmonization = "absoluteChanges")
  )
  expect_true(any(grepl("Check failed.*absorb their own wood harvest area", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional wood harvest area is parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global wood harvest area is parallel", conditions$conditions)))
  x <- result$x

  # regone.1 sits exactly on its cap, its region partner stays below its cap,
  # regtwo absorbs the excess into its headroom, nothing exceeds the land
  aggOut <- toolAggregateWoodHarvest(x[, c(2025, 2030), ])
  maxHarvest <- toolMaxHarvestPerYear(land, split = FALSE)[, c(2025, 2030), ]
  expect_equal(as.vector(aggOut["regone.1", 2025, "secdf"]), 10)
  expect_true(max(aggOut - maxHarvest) <= 10^-10)
  expect_false(any(grepl("Check failed.*wood harvest area is smaller than land", conditions$conditions)))

  # the global total keeps the exact parallelism: regone expects 26, regtwo 4
  regionTotal <- dimSums(dimSums(aggOut, 1.2), 3)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2025, ]), 30, tolerance = 10^-8)
  expect_true(as.vector(regionTotal["regone", 2025, ]) < 26)
})

test_that("absoluteChanges respects the physical cap and reports the unhostable excess", {
  # regone grows its secondary harvest to 9 per cell while both cells only hold
  # 5 Mha of secondary land: not even the global pool has headroom, so the
  # output stays at the cap and the dropped excess is reported as global excess
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 5, secdn = 1, pltns = 1),
    regone.2 = list(primf = 20, primn = 5, secdf = 5, secdn = 1, pltns = 1)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 0.5, secnf = 0.5)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 9, 9), pltns = 0.5, secnf = 0.5),
                             regone.2 = anchor))

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandHarmonized = land,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )
  conditions <- captureConditions(
    result <- calcWoodHarvestAreaHarmonized(input = "magpie", target = "luh3",
                                            harmonizationPeriod = c(2020, 2020),
                                            harmonization = "absoluteChanges")
  )
  # the 4 Mha yr-1 excess of regone.1 is reported for both future years
  expect_true(any(grepl("Check failed.*Harvest area fits into available land", conditions$conditions)))
  expect_true(any(grepl("global excess harvest area", conditions$conditions)))
  expect_true(any(grepl("2025: 4, 2030: 4", conditions$conditions)))
  x <- result$x

  # the output never exceeds the physical cap and stays non-negative
  aggOut <- toolAggregateWoodHarvest(x[, -1, ])
  maxHarvest <- toolMaxHarvestPerYear(land, split = FALSE)
  expect_true(max(aggOut[, c(2025, 2030), ] - maxHarvest[, c(2025, 2030), ]) <= 10^-10)
  expect_true(min(x) >= 0)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 5)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 2)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.pltns"]), 0.5)
})

test_that("fade harmonization keeps its existing behavior", {
  # constant harvest area everywhere: the faded harvest equals the input after
  # the harmonization period, all structural checks stay green, nothing shifts
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 50, primn = 10, secdf = 20, secdn = 10, pltns = 10),
    regone.2 = list(primf = 50, primn = 10, secdf = 20, secdn = 10, pltns = 10),
    regtwo.1 = list(primf = 50, primn = 10, secdf = 20, secdn = 10, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020, 2025), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandHarmonized = land,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )
  conditions <- captureConditions(
    result <- calcWoodHarvestAreaHarmonized(input = "magpie", target = "luh3",
                                            harmonizationPeriod = c(2020, 2030),
                                            harmonization = "fade")
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  expect_true(any(grepl("Harvest area fits into available land", conditions$conditions)))
  x <- result$x

  # after the harmonization period the harvest equals the input
  expect_equal(as.vector(x["regone.1", 2030, "wood_harvest_area.secyf"]), 2)
  expect_equal(as.vector(x["regone.1", 2030, "wood_harvest_area.pltns"]), 1)
  expect_equal(as.vector(x["regone.1", 2030, "wood_harvest_area.secnf"]), 1)
  expect_equal(as.vector(x["regone.1", 2030, "wood_harvest_area.primf"]), 0)
  # the history is untouched
  expect_equal(as.vector(x[, 2010, harvestItems]), as.vector(xTarget[, 2010, harvestItems]))
  expect_true(min(x) >= 0)
})
