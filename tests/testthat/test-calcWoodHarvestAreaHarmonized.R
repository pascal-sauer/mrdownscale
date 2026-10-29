harvestItems <- paste0("wood_harvest_area.", c("primf", "secyf", "secmf", "pltns", "primn", "secnf"))
landItems <- c("primf", "secdf", "pltns", "primn", "secdn")

makeHarvest <- function(years, clusterValues) {
  x <- new.magpie(names(clusterValues), years = years, names = harvestItems, fill = 0)
  getSets(x) <- c("region", "id", "year", "category", "data")
  for (cluster in names(clusterValues)) {
    for (item in names(clusterValues[[cluster]])) {
      x[cluster, , paste0("wood_harvest_area.", item)] <- clusterValues[[cluster]][[item]]
    }
  }
  return(x)
}

makeLand <- function(years, clusterValues) {
  x <- new.magpie(names(clusterValues), years = years, names = landItems, fill = 0)
  getSets(x) <- c("region", "id", "year", "land")
  for (cluster in names(clusterValues)) {
    for (item in names(clusterValues[[cluster]])) {
      x[cluster, , item] <- clusterValues[[cluster]][[item]]
    }
  }
  return(x)
}

test_that("absoluteChanges keeps regional and global harvest area parallel to the input", {
  # regone grows its secondary harvest in the input while its harmonized land
  # stays flat, regtwo is inert: the region and global harvest area totals must
  # equal the target total in the harmonization year plus the regional input
  # delta, the prim harvest stays land-derived (1 Mha yr-1 on regone.1) and the
  # secondary harvest is allocated onto the clusters proportional to their
  # harvestable land
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
  # (8 per region with two clusters, 4 in regtwo) plus the input delta
  regionTotal <- dimSums(dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2), 3)
  expect_equal(as.vector(regionTotal["regone", 2025, ]), 8 + (14 - 8))
  expect_equal(as.vector(regionTotal["regone", 2030, ]), 8 + (20 - 8))
  expect_equal(as.vector(regionTotal["regtwo", 2025, ]), 4)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2025, ]), 14 + 4)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2030, ]), 20 + 4)

  # the prim harvest stays exactly land-derived, the rest of the forest cover
  # target (11 in 2025) is allocated onto clusters proportional to their
  # harvestable forest land (60:50), the untouched other land cover stays on
  # its anchor of 2
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.primf"]), 1, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 11 * 50 / 110, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 11 * 40 / 110, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.pltns"]), 11 * 10 / 110, tolerance = 10^-8)
  # inert regtwo keeps its cover totals, allocated onto its harvestable land
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secyf"]), 3 * 30 / 40, tolerance = 10^-8)
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.pltns"]), 3 * 10 / 40, tolerance = 10^-8)
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secnf"]), 1)
  # the target history is untouched
  expect_equal(as.vector(x[, 2020, harvestItems]), as.vector(xTarget[, 2020, harvestItems]))
})

test_that("absoluteChanges shares a cover harvest within a region proportional to harvestable land", {
  # regone grows its forest harvest on regone.1 while its harmonized forest
  # land stays flat: the regional cover target is allocated onto the clusters
  # proportional to their harvestable forest land (secdf 60:20 -> 3:1), not to
  # the raw input shares, so the cluster with three times the land gets three
  # times the harvest
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 10, primn = 5, secdf = 60, secdn = 10, pltns = 0),
    regone.2 = list(primf = 10, primn = 5, secdf = 20, secdn = 5, pltns = 0)
  ))

  xTarget <- makeHarvest(c(2010, 2020), list(
    regone.1 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 0, secnf = 1),
    regone.2 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 0, secnf = 1)
  ))
  xInput <- makeHarvest(c(2020, 2025, 2030), list(
    regone.1 = list(secyf = c(2, 10, 10), pltns = 0, secnf = 1),
    regone.2 = list(secyf = 2, pltns = 0, secnf = 1)
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

  # the forest target 4 + (12 - 4) = 12 is shared 9:3 by the 3:1 caps
  aggOut <- toolAggregateWoodHarvest(x[, c(2025, 2030), ])
  expect_equal(as.vector(aggOut["regone.1", 2025, "secdf"]), 9, tolerance = 10^-8)
  expect_equal(as.vector(aggOut["regone.2", 2025, "secdf"]), 3, tolerance = 10^-8)
  # the other land cover target 2 follows its harvestable land (10:5)
  expect_equal(as.vector(aggOut["regone.1", 2025, "secdn"]), 4 / 3, tolerance = 10^-8)
  expect_equal(as.vector(aggOut["regone.2", 2025, "secdn"]), 2 / 3, tolerance = 10^-8)
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges keeps regional forest and other land cover totals parallel to the input", {
  # regone grows its forest secondary harvest while its other land harvest
  # collapses below the target anchor: the target is broken down per cover
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

  # the cover targets, rebuilt from the mocks exactly the way the
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
    target <- setYears(dimSums(targetRegion[, 2020, items], 3), NULL) +
      (dimSums(inputRegion[, c(2025, 2030), items], 3) -
         setYears(dimSums(inputRegion[, 2020, items], 3), NULL))
    target[target < 0] <- 0
    actual <- dimSums(regionTotal[, , items], 3)
    # in 2025: regone forest 6 + 8 = 14, other 1.2 - 1 = 0.2, regtwo 3 and 1
    expect_equal(as.data.frame(actual), as.data.frame(target), tolerance = 10^-8)
  }
  # the global cover totals stay parallel as well
  expect_equal(as.vector(dimSums(dimSums(regionTotal[, , covers$forest], 3), 1)[, 2025, ]), 17,
               tolerance = 10^-8)
  expect_equal(as.vector(dimSums(dimSums(regionTotal[, , covers$other], 3), 1)[, 2025, ]), 1.2,
               tolerance = 10^-8)
})

test_that("absoluteChanges shares a collapsed cover target over the region's harvestable land", {
  # regone.1's raw absolute change turns negative and is floored to zero, so
  # the region forest target falls short of the summed cluster changes: the
  # target 3.5 is allocated onto the region's harvestable forest land
  # (60:50), giving the collapsing cluster its proportional share of 1.59 instead
  # of the zero the raw clamping would have left it with, the other land cover
  # is untouched, inert regtwo keeps its target of 3 on its own land
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regone.2 = list(primf = 20, primn = 5, secdf = 40, secdn = 10, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 15, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  # regone.1's target is only 1 but its input loses 1.5 from 2 to 0.5, so its
  # raw absolute change of -0.5 turns negative and gets floored to 0
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

  # the forest target 3.5 is shared proportional to the harvestable
  # forest land (50:40 secdf), its region partner no longer compensates alone
  expect_equal(as.vector(x["regone.1", 2025, "wood_harvest_area.secyf"]), 3.5 * 50 / 110, tolerance = 10^-8)
  expect_equal(as.vector(x["regone.2", 2025, "wood_harvest_area.secyf"]), 3.5 * 40 / 110, tolerance = 10^-8)
  expect_equal(as.vector(x["regtwo.1", 2025, "wood_harvest_area.secyf"]), 3 * 30 / 40, tolerance = 10^-8)

  # the cover totals stay parallel: forest 5 + (4.5 - 6) = 3.5, other 2
  regionTotal <- dimSums(dimSums(toolAggregateWoodHarvest(x[, c(2025, 2030), ]), 1.2), 3)
  expect_equal(as.vector(regionTotal["regone", 2025, ]), 5.5, tolerance = 10^-8)
  expect_equal(as.vector(regionTotal["regtwo", 2025, ]), 4)
})

test_that("absoluteChanges flows a region target above its cap into the other cover", {
  # regone targets 34 aggregated forest harvest in 2025 but its forest land
  # only allows 25: both clusters are allocated exactly on their caps, the
  # remainder of 9 fills the headroom of regone's other land cover (8) and,
  # over the global margin, regtwo's other land headroom (1): both regional
  # cover checks warn, the caps hold by construction and the global total
  # stays exactly parallel while no excess is reported
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
  expect_true(any(grepl("Check failed.*Regional forest wood harvest area is parallel", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional other wood harvest area is parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global wood harvest area is parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Harvest area fits into available land", conditions$conditions)))
  x <- result$x

  # both regone forest clusters sit exactly on their caps, the transferred other
  # land harvest fills their other land clusters exactly to their caps as well
  aggOut <- toolAggregateWoodHarvest(x[, c(2025, 2030), ])
  maxHarvest <- toolMaxHarvestPerYear(land, split = FALSE)[, c(2025, 2030), ]
  expect_equal(as.vector(aggOut["regone.1", 2025, c("secdf", "pltns")]), c(10, 5))
  expect_equal(as.vector(aggOut["regone.2", 2025, c("secdf", "pltns")]), c(5, 5))
  expect_equal(as.vector(aggOut["regone.1", 2025, "secdn"]), 5)
  expect_equal(as.vector(aggOut["regone.2", 2025, "secdn"]), 5)
  # regtwo's forest target 3 lands proportional to its forest land (30:5)
  expect_equal(as.vector(aggOut["regtwo.1", 2025, c("secdf", "pltns")]),
               c(3 * 30 / 35, 3 * 5 / 35), tolerance = 10^-8)
  # its other land cluster holds its target 1 plus the 1 received globally
  expect_equal(as.vector(aggOut["regtwo.1", 2025, "secdn"]), 2, tolerance = 10^-8)
  expect_true(max(aggOut - maxHarvest) <= 10^-10)
  expect_false(any(grepl("Check failed.*wood harvest area is smaller than land", conditions$conditions)))

  # the global total keeps the exact parallelism: regone targets 36, regtwo 4
  regionTotal <- dimSums(dimSums(aggOut, 1.2), 3)
  expect_equal(as.vector(dimSums(regionTotal, 1)[, 2025, ]), 40, tolerance = 10^-8)
  expect_equal(as.vector(regionTotal["regone", 2025, ]), 35, tolerance = 10^-8)
})

test_that("absoluteChanges respects the physical cap and reports the unhostable excess", {
  # regone grows its forest harvest to 20 on regone.1 while the forest land of
  # both clusters together only allows 12 and the 1 Mha yr-1 of other land
  # headroom takes just a part of the 11 excess: what even the other cover
  # cannot host is dropped at the cap and reported as global excess of 10
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
  # the 10 Mha yr-1 excess that not even the other cover can host is reported
  # for both future years, the warn-only parallel checks report the shortfall
  expect_true(any(grepl("Check failed.*Harvest area fits into available land", conditions$conditions)))
  expect_true(any(grepl("global excess harvest area", conditions$conditions)))
  expect_true(any(grepl("2025: 10, 2030: 10", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional forest wood harvest area is parallel", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional other wood harvest area is parallel", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Global wood harvest area is parallel", conditions$conditions)))
  x <- result$x

  # both clusters are allocated exactly on their caps, the hostable other land
  # target fills its clusters to their caps, nothing exceeds the land
  aggOut <- toolAggregateWoodHarvest(x[, -1, ])
  maxHarvest <- toolMaxHarvestPerYear(land, split = FALSE)
  expect_true(max(aggOut[, c(2025, 2030), ] - maxHarvest[, c(2025, 2030), ]) <= 10^-10)
  expect_true(min(x) >= 0)
  expect_equal(as.vector(aggOut["regone.1", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(aggOut["regone.2", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(aggOut["regone.1", 2025, "secdn"]), 1)
  expect_equal(as.vector(aggOut["regone.2", 2025, "secdn"]), 1)
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

test_that("absoluteChanges allocates cover harvest onto clusters without own harvest of that cover", {
  # regone.2 has no forest harvest at all but forest land, while regone.1's
  # forest target exceeds its own cap: the allocation lands on both clusters
  # proportional to their forest caps, regone.1 exactly on its 6 Mha cap and
  # regone.2 filling its 15 Mha with nothing to start from
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

  # the forest cover stays region parallel: regone's 21 Mha yr-1 target
  # fits exactly into the region's 21 Mha of forest land, regone.1 sitting on
  # its 6 Mha cap and regone.2 filling its 15 Mha with nothing to start from
  aggOut <- toolAggregateWoodHarvest(x[, c(2025, 2030), ])
  expect_equal(as.vector(dimSums(dimSums(aggOut, 1.2)[, , c("primf", "secdf", "pltns")], 3)["regone", 2025, ]),
               21, tolerance = 10^-8)
  expect_equal(as.vector(aggOut["regone.1", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(aggOut["regone.2", 2025, c("secdf", "pltns")]), c(10, 5))
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges keeps the global cover totals parallel when a regional target is negative", {
  # regone's other land input collapses far below its target anchor, so its
  # un-floored other land target turns negative while regtwo grows: the
  # floored regional target must not lift the global cover total above the
  # raw global target, regtwo is scaled down to keep it exactly parallel
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 15, pltns = 10)
  ))

  # the forest cover stays flat at its anchor in both regions, its targets
  # never hit a floor and must stay exactly as they are today
  xTarget <- makeHarvest(c(2010, 2020), list(
    regone.1 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4),
    regtwo.1 = list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4)
  ))
  # the un-floored other land targets are regone -5 and regtwo 3 in 2025,
  # regone -5 and regtwo 12 in 2030, giving the raw global targets 3 and 7
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

  # the global other land total equals the un-floored raw global target
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(otherTotal, 1)[, 2030, ]), 7, tolerance = 10^-9)

  # the loss-making region gets nothing of the cover
  expect_equal(as.vector(otherTotal["regone", 2025, ]), 0)
  expect_equal(as.vector(otherTotal["regone", 2030, ]), 0)

  # the growing region carries the whole global target, scaled below its own
  expect_equal(as.vector(otherTotal["regtwo", 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(otherTotal["regtwo", 2030, ]), 7, tolerance = 10^-9)

  # the forest cover with its unharmed targets stays exactly the same
  forestTotal <- dimSums(regionTotal[, , c("primf", "secdf", "pltns")], 3)
  expect_equal(as.vector(forestTotal["regone", 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(forestTotal["regtwo", 2030, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2025, ]), 6, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2030, ]), 6, tolerance = 10^-9)

  expect_true(min(x) >= 0)
})

test_that("absoluteChanges charges a negative global cover target to the other cover", {
  # both regions' other land input collapses so that the global other land
  # target turns negative while the global forest target outweighs
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
  # the un-floored cover targets are forest 3 + 3 = 6 in both regions and
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
  # un-floored global target 6 - 3, not its own global target 6
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal, 1)[, 2030, ]), 3, tolerance = 10^-9)

  # the global total over both covers equals the un-floored target 3
  expect_equal(as.vector(dimSums(forestTotal + otherTotal, 1)[, 2025, ]), 3, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forestTotal + otherTotal, 1)[, 2030, ]), 3, tolerance = 10^-9)

  # both positive forest regions are scaled by the balancing factor 3 / 6
  expect_equal(as.vector(forestTotal["regone", 2025, ]), 1.5, tolerance = 10^-9)
  expect_equal(as.vector(forestTotal["regtwo", 2025, ]), 1.5, tolerance = 10^-9)
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges cross-cover balancing is a no-op without a negative global cover target", {
  # regone's other land target collapses below zero but regtwo's growth
  # keeps the global other land target positive: the balancing factor is
  # exactly 1 and both cover global totals stay parallel as before
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 5, secdf = 50, secdn = 20, pltns = 10),
    regtwo.1 = list(primf = 20, primn = 5, secdf = 30, secdn = 20, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 4)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regtwo.1 = anchor))
  # the un-floored cover targets are forest 3 + 3 = 6 and other -5 in
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

  # both cover global totals equal their own un-floored global targets
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
