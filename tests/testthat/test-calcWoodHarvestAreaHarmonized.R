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
  # harvest already over-achieves scales down both regone forest cells, the
  # untouched other land cover stays on its anchor
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.primf"]), 1, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 8 * 11 / 12, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 2 * 11 / 12, tolerance = 10^-8)
  # inert regtwo keeps its anchor, the target history is untouched
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secyf"]), 2)
  expect_equal(as.vector(x[, 2020, harvestItems]), as.vector(xTarget[, 2020, harvestItems]))
})

test_that("absoluteChanges keeps regional forest and other land cover totals parallel to the input", {
  # regone grows its forest secondary harvest while its other land harvest
  # collapses below the target anchor: the expectation is broken down per cover
  # (forest = primf, secdf, pltns, other = primn, secdn) and both cover totals
  # must stay parallel within the region, which the pooled rescaling cannot
  # achieve because there the forest growth absorbs the other land loss
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regone.2 = list(primf = 20, primn = 5, secdf = 40, secdn = 10, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 15, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  # regone.1 grows its forest harvest from 2 to 10 (2025) and 16 (2030) while
  # its other land harvest collapses from 1 to 0, far below its small target
  # anchor of 0.2, so the raw absolute change clamps to zero inside one cover
  xTarget["regone.1", , "wood_harvest_area.secnf"] <- 0.2
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 10, 16), pltns = 1, secnf = c(1, 0, 0)),
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

  # the expected cover totals, rebuilt from the mocks exactly the way the
  # module promises them: floored target cover total in 2020 plus the regional
  # input cover delta, the flat land leaves the prim harvest at zero and the
  # ample secondary land never binds
  covers <- list(forest = c("primf", "secdf", "pltns"), other = c("primn", "secdn"))
  regionTotal <- dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2)
  stopifnot(all(unlist(covers) %in% getItems(regionTotal, 3)))
  inputRegion <- dimSums(toolAggregateWoodHarvest(xInput), 1.2)
  targetRegion <- dimSums(toolAggregateWoodHarvest(xTarget), 1.2)
  for (cover in names(covers)) {
    items <- covers[[cover]]
    expected <- setYears(dimSums(targetRegion[, 2020, items], 3), NULL) +
      (dimSums(inputRegion[, c(2025, 2030), items], 3) -
         setYears(dimSums(inputRegion[, 2020, items], 3), NULL))
    expected[expected < 0] <- 0
    actual <- dimSums(regionTotal[, , items], 3)
    # in 2025: regone forest 6 + 8 = 14, other 1.2 - 1 = 0.2, regtwo 3 and 1
    expect_equal(as.data.frame(actual), as.data.frame(expected), tolerance = 10^-8)
  }
  # the global cover totals stay parallel as well
  expect_equal(as.vector(dimSums(dimSums(regionTotal[, , covers$forest], 3), 1)[, 2025, ]), 17,
               tolerance = 10^-8)
  expect_equal(as.vector(dimSums(dimSums(regionTotal[, , covers$other], 3), 1)[, 2025, ]), 1.2,
               tolerance = 10^-8)
})

test_that("absoluteChanges compensates clamped losses within the same region", {
  # regone.1's raw absolute change turns negative and is clamped to zero, so
  # its region forest total falls short of the expectation and must be
  # compensated by the other cells of the same forest cover: the forest
  # expectation 3.5 over an actual 4 scales regone.2 from 2 to 1.75, the other
  # land cover is untouched, inert regtwo keeps its anchor of 2
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
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 2 * 3.5 / 4, tolerance = 10^-8)
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secyf"]), 2)

  # the cover totals stay parallel: forest 5 + (4.5 - 6) = 3.5, other 2
  regionTotal <- dimSums(dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2), 3)
  expect_equal(as.vector(regionTotal["regone", 2025, ]), 5.5, tolerance = 10^-8)
  expect_equal(as.vector(regionTotal["regtwo", 2025, ]), 4)
})

test_that("absoluteChanges falls back to the global pool if a region cannot host the harvest", {
  # regone expects 34 aggregated forest harvest in 2025 but its forest land
  # only allows 25, so the excess of regone.1 goes over the global pool into
  # regtwo's forest headroom: the warn check fires, the forest cover is
  # off-parallel, the untouched other land cover stays parallel, but the global
  # total stays parallel while the global capacity suffices
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 10, secdn = 5, pltns = 5),
    regone.2 = list(primf = 20, primn = 5, secdf = 5, secdn = 5, pltns = 5),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 5, pltns = 5)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 30, 30), pltns = 1, secnf = 1),
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
  expect_true(any(grepl("Check failed.*absorb their own.*wood harvest area", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional forest wood harvest area is parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Regional other wood harvest area is parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global wood harvest area is parallel", conditions$conditions)))
  x <- result$x

  # regone.1 sits exactly on its cap, its region partner stays below its cap,
  # regtwo absorbs the excess into its headroom, nothing exceeds the land
  aggOut <- toolAggregateWoodHarvest(x[, c(2025, 2030), ])
  maxHarvest <- toolMaxHarvestPerYear(land, split = FALSE)[, c(2025, 2030), ]
  expect_equal(as.vector(aggOut["regone.1", 2025, "secdf"]), 10)
  expect_true(max(aggOut - maxHarvest) <= 10^-10)
  expect_false(any(grepl("Check failed.*wood harvest area is smaller than land", conditions$conditions)))

  # the global total keeps the exact parallelism: regone expects 36, regtwo 4
  regionTotal <- dimSums(dimSums(aggOut, 1.2), 3)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2025, ]), 40, tolerance = 10^-8)
  expect_true(as.vector(regionTotal["regone", 2025, ]) < 36)
})

test_that("absoluteChanges respects the physical cap and reports the unhostable excess", {
  # regone grows its forest harvest to 20 on regone.1 while the forest land of
  # both cells together only allows 12 and the 1 Mha yr-1 of other land
  # headroom takes just a part of the excess: what even the other cover cannot
  # host is dropped at the cap and reported as global excess
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 5, secdn = 1, pltns = 1),
    regone.2 = list(primf = 20, primn = 5, secdf = 5, secdn = 1, pltns = 1)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 0.5, secnf = 0.5)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 20, 20), pltns = 0.5, secnf = 0.5),
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
  # the 13.5 Mha yr-1 excess of regone.1 is reported for both future years
  expect_true(any(grepl("Check failed.*Harvest area fits into available land", conditions$conditions)))
  expect_true(any(grepl("global excess harvest area", conditions$conditions)))
  expect_true(any(grepl("2025: 13.5, 2030: 13.5", conditions$conditions)))
  x <- result$x

  # the output never exceeds the physical cap and stays non-negative
  aggOut <- toolAggregateWoodHarvest(x[, -1, ])
  maxHarvest <- toolMaxHarvestPerYear(land, split = FALSE)
  expect_true(max(aggOut[, c(2025, 2030), ] - maxHarvest[, c(2025, 2030), ]) <= 10^-10)
  expect_true(min(x) >= 0)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 5)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 2)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.pltns"]), 1)
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

test_that("absoluteChanges can move cover harvest into cells without own harvest of that cover", {
  # regone.2 has no forest harvest at all but forest land, while regone.1's
  # forest excess exceeds its own cap: the region absorbs the excess into
  # regone.2's forest headroom, so the cover total must be filled onto an item
  # split of a cell that held nothing of that cover before
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 5, secdn = 1, pltns = 1),
    regone.2 = list(primf = 20, primn = 5, secdf = 10, secdn = 1, pltns = 5)
  ))

  xTarget <- makeHarvest(c(2010, 2020), list(
    regone.1 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 0),
    regone.2 = list(primf = 0, primn = 0, secyf = 0, secmf = 0, pltns = 0, secnf = 0)
  ))
  xInput <- makeHarvest(c(2020, 2025, 2030), list(
    regone.1 = list(secyf = c(2, 20, 20), pltns = 1, secnf = 0),
    regone.2 = list(secnf = 0)
  ))

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

  # the forest cover stays region parallel: regone's 21 Mha yr-1 expectation
  # fits exactly into the region's 21 Mha of forest land, regone.1 sitting on
  # its 6 Mha cap and regone.2 filling its 15 Mha with nothing to start from
  aggOut <- toolAggregateWoodHarvest(x[, c(2025, 2030), ])
  expect_equal(as.vector(dimSums(dimSums(aggOut, 1.2)[, , c("primf", "secdf", "pltns")], 3)["regone", 2025, ]),
               21, tolerance = 10^-8)
  expect_equal(as.vector(aggOut["regone.1", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(aggOut["regone.2", 2025, c("secdf", "pltns")]), c(10, 5))
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges keeps the global cover totals parallel when a regional expectation is negative", {
  # regone's other land input collapses far below its target anchor, so its
  # un-floored other land expectation turns negative while regtwo grows: the
  # floored regional expectation must not lift the global cover total above the
  # raw global expectation, regtwo is scaled down to keep it exactly parallel
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 15, pltns = 10)
  ))

  # the forest cover stays flat at its anchor in both regions, its expectations
  # never hit a floor and must stay exactly as they are today
  xTarget <- makeHarvest(c(2010, 2020), list(
    regone.1 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4),
    regtwo.1 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4)
  ))
  # the un-floored other land expectations are regone -5 and regtwo 3 in 2025,
  # regone -5 and regtwo 12 in 2030, giving the raw global expectations 3 and 7
  xInput <- makeHarvest(c(2020, 2025, 2030), list(
    regone.1 = list(secyf = 2, pltns = 1, secnf = c(10, 1, 1)),
    regtwo.1 = list(secyf = 2, pltns = 1, secnf = c(4, 8, 12))
  ))

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

  regionTotal <- dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2)
  otherTotal <- dimSums(regionTotal[, , c("primn", "secdn")], 3)

  # the global other land total equals the un-floored raw global expectation
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2030, ]), 7, tolerance = 10^-9)

  # the loss-making region gets nothing of the cover
  expect_equal(as.vector(otherTotal["regone", 2025, ]), 0)
  expect_equal(as.vector(otherTotal["regone", 2030, ]), 0)

  # the growing region carries the whole global expectation, scaled below its own
  expect_equal(as.vector(otherTotal["regtwo", 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(otherTotal["regtwo", 2030, ]), 7, tolerance = 10^-9)

  # the forest cover with its unharmed expectations stays exactly the same
  forestTotal <- dimSums(regionTotal[, , c("primf", "secdf", "pltns")], 3)
  expect_equal(as.vector(forestTotal["regone", 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(forestTotal["regtwo", 2030, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2025, ]), 6, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2030, ]), 6, tolerance = 10^-9)

  expect_true(min(x) >= 0)
})

test_that("absoluteChanges charges a negative global cover expectation to the other cover", {
  # both regions' other land input collapses so that the global other land
  # expectation turns negative while the global forest expectation outweighs
  # the loss: the floored other land total must be charged to the forest cover
  # so the global total over all covers stays parallel to the input and the
  # positive forest regions are scaled by the balancing factor 3 / 6
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 15, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regtwo.1 = anchor))
  # the un-floored cover expectations are forest 3 + 3 = 6 in both regions and
  # other -5 in regone plus 2 in regtwo, giving the global totals 6 and -3
  xInput <- makeHarvest(c(2020, 2025, 2030), list(
    regone.1 = list(secyf = 2, pltns = 1, secnf = c(10, 1, 1)),
    regtwo.1 = list(secyf = 2, pltns = 1, secnf = c(4, 2, 2))
  ))

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

  regionTotal <- dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2)
  forestTotal <- dimSums(regionTotal[, , c("primf", "secdf", "pltns")], 3)
  otherTotal <- dimSums(regionTotal[, , c("primn", "secdn")], 3)

  # the floored other land cover gives up its whole global total
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2025, ]), 0)
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2030, ]), 0)

  # the forest cover absorbs the charged loss: its global total is the summed
  # un-floored global expectation 6 - 3, not its own global expectation 6
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2030, ]), 3, tolerance = 10^-9)

  # the global total over both covers equals the un-floored expectation 3
  expect_equal(as.vector(dimSums(forestTotal + otherTotal, 1)[, 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal + otherTotal, 1)[, 2030, ]), 3, tolerance = 10^-9)

  # both positive forest regions are scaled by the balancing factor 3 / 6
  expect_equal(as.vector(forestTotal["regone", 2025, ]), 1.5, tolerance = 10^-9)
  expect_equal(as.vector(forestTotal["regtwo", 2025, ]), 1.5, tolerance = 10^-9)
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges cross-cover balancing is a no-op without a negative global cover expectation", {
  # regone's other land expectation collapses below zero but regtwo's growth
  # keeps the global other land expectation positive: the balancing factor is
  # exactly 1 and both cover global totals stay parallel as before
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 20, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regtwo.1 = anchor))
  # the un-floored cover expectations are forest 3 + 3 = 6 and other -5 in
  # regone plus 20 in regtwo, giving the positive global totals 6 and 15
  xInput <- makeHarvest(c(2020, 2025, 2030), list(
    regone.1 = list(secyf = 2, pltns = 1, secnf = c(10, 1, 1)),
    regtwo.1 = list(secyf = 2, pltns = 1, secnf = c(4, 20, 20))
  ))

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

  regionTotal <- dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2)
  forestTotal <- dimSums(regionTotal[, , c("primf", "secdf", "pltns")], 3)
  otherTotal <- dimSums(regionTotal[, , c("primn", "secdn")], 3)

  # both cover global totals equal their own un-floored global expectations
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2025, ]), 6, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2030, ]), 6, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2025, ]), 15, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2030, ]), 15, tolerance = 10^-9)

  # the loss-making region stays at zero, the growing region keeps its floored
  # and scaled total, its region partner stays parallel
  expect_equal(as.vector(otherTotal["regone", 2025, ]), 0)
  expect_equal(as.vector(otherTotal["regtwo", 2025, ]), 15, tolerance = 10^-9)
  expect_equal(as.vector(forestTotal["regone", 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(forestTotal["regtwo", 2025, ]), 3, tolerance = 10^-9)
  expect_true(min(x) >= 0)
})
