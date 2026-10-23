nonlandItems <- c("bioh.primf", "bioh.secmf", "wood_harvest_area.primf", "wood_harvest_area.secmf",
                  "harvest_weight_type.roundwood", "harvest_weight_type.fuelwood", "fertilizer.c3ann")

test_that("calcNonlandHarmonized works with harmonization = absoluteChanges", {
  harmonizationYear <- 2020

  # real MAgPIE/LUH3 source data is not available in CI, so the source data
  # functions are mocked with synthetic data
  # for reg.one the kg C per Mha of primf drops in the input from 2000 to 500,
  # the raw absolute change is negative, so the harmonizer clamps it to zero
  # and it is later floored to the smallest positive rate
  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025, 2030), names = nonlandItems, fill = 0)
  xInput["reg.one", 2020, ] <- c(2000, 1500, 1, 1, 3, 1, 400)
  xInput["reg.one", 2025, ] <- c(500, 1600, 1, 1, 2, 2, 500)
  # in 2030 reg.two has cropland again, so no fertilizer is lost that year
  xInput["reg.one", 2030, ] <- c(500, 1600, 1, 1, 2, 2, 450)
  xInput["reg.two", 2020, ] <- c(1000, 1000, 1, 1, 1, 1, 500)
  xInput["reg.two", 2025, ] <- c(1200, 1000, 1, 1, 1, 1, 600)
  xInput["reg.two", 2030, ] <- c(1200, 1000, 1, 1, 1, 1, 550)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget["reg.one", 2010, ] <- c(1000, 1500, 1, 1, 2, 2, 500)
  xTarget["reg.one", 2020, ] <- c(1000, 1500, 1, 1, 2, 2, 500)
  xTarget["reg.two", 2010, ] <- c(1000, 1000, 1, 1, 1, 1, 400)
  xTarget["reg.two", 2020, ] <- c(1000, 1000, 1, 1, 1, 1, 400)
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  # fertilizer in the target data is in kg ha-1 yr-1, with 100 Mha of cropland
  # it corresponds to 50 Tg in reg.one and 40 Tg in reg.two
  landTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020),
                           names = c("c3ann", "pastr"), fill = 100)
  landHarmonized <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020, 2025, 2030),
                               names = c("c3ann", "pastr"), fill = 100)
  # in 2025 reg.one raw fertilizer of 150 Tg divided by 100 Mha gives 1500 kg
  # ha-1 yr-1, reg.two has no cropland left, turning its raw fertilizer of
  # 140 Tg into an infinite rate, which must be set to zero
  landHarmonized["reg.two", 2025, "c3ann"] <- 0

  # wood harvest area stays at 1 Mha yr-1 everywhere
  harvestArea <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020, 2025, 2030),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  # temporarily replaces calcOutput in the package namespace, so the internal
  # source data calls return the synthetic data above and anything else errors
  # for equal harmonization years the Extrapolated targets are the LowRes
  # targets returned unchanged (see test-calcNonlandTargetExtrapolated.R)
  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  # the rescaled fertilizer rate of 2900 kg ha-1 yr-1 exceeds the plausibility
  # threshold of toolCheckFertilizer, which reports it as a failed check
  expect_no_warning(
    expect_message(
      result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                      harmonizationPeriod = c(harmonizationYear, harmonizationYear),
                                      harmonization = "absoluteChanges"),
      "Fertilizer application"
    )
  )
  x <- result$x

  expect_false(is.null(attr(x, "geometry")))
  expect_false(is.null(attr(x, "crs")))
  expect_true(setequal(getItems(x, 3), nonlandItems))
  expect_equal(getYears(x, as.integer = TRUE), c(2010, 2020, 2025, 2030))

  # before and in the harmonization year target data is returned unchanged
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)

  # infinite fertilizer rates are replaced with zero, the resulting global loss
  # is compensated by rescaling the rates of the cells with remaining cropland,
  # so the global totals stay parallel to the input: target total in 2020
  # (50 + 40 Tg) plus global input delta -> 290 Tg in 2025, reached by scaling
  # reg.one's 150 Tg by 290/150 to 2900 kg ha-1 yr-1 on its 100 Mha of cropland
  expect_equal(as.vector(x["reg.one", 2025, "fertilizer"]), 2900)
  expect_equal(as.vector(x["reg.two", 2025, "fertilizer"]), 0)
  # in 2030 no fertilizer is lost, so the absolute changes are reproduced exactly
  expect_equal(as.vector(x["reg.one", 2030, "fertilizer"]), 1000)
  expect_equal(as.vector(x["reg.two", 2030, "fertilizer"]), 900)
  expect_true(all(is.finite(x)))

  fertilizerGlobalTg <- dimSums(toolFertilizerTg(x[, c(2025, 2030), "fertilizer"], landHarmonized), 1)
  expect_equal(as.vector(fertilizerGlobalTg[, 2025, ]), (50 + 40) + (1100 - 900))
  expect_equal(as.vector(fertilizerGlobalTg[, 2030, ]), (50 + 40) + (1000 - 900))

  # the clamped kg C per Mha rate of 0 for reg.one primf is floored to the
  # smallest positive rate (1000), bioh is recalculated from it and renormalized
  # so that the total bioh of 1600 kg C yr-1 of the harmonizer output is kept:
  # primf gets 1600 * 1000 / 2600, secmf 1600 * 1600 / 2600
  expect_equal(as.vector(x["reg.one", 2025, c("bioh.primf", "bioh.secmf")]),
               c(1600 * 1000 / 2600, 1600 * 1600 / 2600))

  # where nothing was clamped, the absolute changes are reproduced exactly
  expect_equal(as.vector(x["reg.two", 2025, c("bioh.primf", "bioh.secmf")]), c(1200, 1000))
  expect_equal(as.vector(x["reg.one", 2025, c("harvest_weight_type.roundwood",
                                              "harvest_weight_type.fuelwood")]), c(1, 3))

  # wood harvest area equals the harmonized harvest area
  expect_equal(as.vector(x["reg.one", 2025, "wood_harvest_area.primf"]), 1)
  expect_true(all(x >= 0))
})

test_that("toolCheckFertilizer also runs for absoluteChanges", {
  # the absolute changes of the input fertilizer (100 to 300 Tg on top of the
  # 5 Tg of the target in the harmonization year) divided by 100 Mha of
  # cropland give 2050 kg ha-1 yr-1 in 2025, which is above the plausibility
  # threshold of 1200 that toolCheckFertilizer reports as a failed check
  xInput <- new.magpie("reg.one", years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, 2020, ] <- c(1000, 1000, 1, 1, 1, 1, 100)
  xInput[, 2025, ] <- c(1000, 1000, 1, 1, 1, 1, 300)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie("reg.one", years = c(2010, 2020), names = nonlandItems, fill = 1)
  xTarget[, , "fertilizer.c3ann"] <- 50
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  land <- new.magpie("reg.one", years = c(2010, 2020, 2025), names = c("c3ann", "pastr"), fill = 100)

  harvestArea <- new.magpie("reg.one", years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = land[, c(2010, 2020), ],
             LandHarmonized = land,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  expect_message(
    calcNonlandHarmonized(input = "magpie", target = "luh3",
                          harmonizationPeriod = c(2020, 2020),
                          harmonization = "absoluteChanges"),
    "Fertilizer application"
  )
})

test_that("absoluteChanges projects input deltas onto the harmonized land", {
  # regone.1 keeps its input fertilizer growth in MAgPIE while its harmonized
  # cropland collapses after 2020, the region's other cell expands: under
  # absoluteChanges the regional tonnage delta must follow the harmonized land,
  # so the collapsing cell keeps only its target anchor plus its tiny area
  # share, and no fertilizer ends up on land that vanished. The tonnages are
  # scaled so the cap below does not interfere with the projection.
  cells <- c("regone.1", "regone.2", "regtwo.1")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput["regone.1", , "fertilizer.c3ann"] <- c(10, 14)
  xInput["regone.2", , "fertilizer.c3ann"] <- c(1, 1)
  xInput["regtwo.1", , "fertilizer.c3ann"] <- c(1, 1)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  # target cropland of 0.01, 20 and 10 Mha anchors the target rate of
  # 100 kg ha-1 yr-1 at 0.001, 2 and 1 Tg
  xTarget <- new.magpie(cells, years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 10)
  landTarget["regone.1", , ] <- 0.01
  landTarget["regone.2", , ] <- 20
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = "c3ann", fill = 10)
  landHarmonized["regone.1", , ] <- 0.01
  landHarmonized["regone.2", , ] <- 20
  landHarmonized["regone.1", 2025, ] <- 0.003

  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  conditions <- captureConditions(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(2020, 2020),
                                    harmonization = "absoluteChanges")
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  x <- result$x

  # the +4 Tg regional delta is distributed proportional to the 2025 harmonized
  # areas of the region (0.003 and 20 Mha), the cell in the region without any
  # delta keeps its anchor rate
  expect_equal(as.vector(x["regone.1", 2025, "fertilizer"]),
               (0.001 + 4 * (0.003 / 20.003)) / 0.003 * (10^9 / 10^6), tolerance = 10^-9)
  expect_equal(as.vector(x["regone.2", 2025, "fertilizer"]),
               (2 + 4 * (20 / 20.003)) / 20 * (10^9 / 10^6), tolerance = 10^-9)
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer"]), 100)
  expect_true(max(x[, 2025, "fertilizer"]) <= 1200)

  # the projection preserves the regional tonnage sums: target plus regional
  # input delta, also globally (which keeps the parallelism check exact)
  regionTg <- dimSums(toolFertilizerTg(x[, 2025, "fertilizer"], landHarmonized), 1.2)
  expect_equal(as.vector(regionTg["regone", , ]), 2.001 + 4)
  expect_equal(as.vector(regionTg["regtwo", , ]), 1)
  expect_equal(as.vector(dimSums(toolFertilizerTg(x[, 2025, "fertilizer"], landHarmonized), c(1, 2))), 3.001 + 4)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
})

test_that("the residual cap pulls collapsed-area cells with anchor tonnage back to 1200", {
  # with a zero regional delta the collapsing cell regone.1 keeps its anchor of
  # 10 Tg on 1e-5 Mha of harmonized cropland, a rate of 1e9 kg ha-1 yr-1: the
  # cap must pull it back to the threshold and redistribute its excess to the
  # headroom of its region partner, leaving the other region untouched
  cells <- c("regone.1", "regone.2", "regtwo.1")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput["regone.1", , "fertilizer.c3ann"] <- c(10, 10)
  xInput["regone.2", , "fertilizer.c3ann"] <- c(10, 10)
  xInput["regtwo.1", , "fertilizer.c3ann"] <- c(1, 1)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget["regone.1", , "fertilizer.c3ann"] <- 1000
  xTarget["regone.2", , "fertilizer.c3ann"] <- 100
  xTarget["regtwo.1", , "fertilizer.c3ann"] <- 10
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  # anchors: 10 Tg in regone.1 (1000 kg ha-1 on 10 Mha) and regone.2
  # (100 kg ha-1 on 100 Mha), 1 Tg in regtwo
  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landTarget["regone.1", , ] <- 10
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = "c3ann", fill = 100)
  landHarmonized["regone.1", c(2010, 2020), ] <- 10
  landHarmonized["regone.1", 2025, ] <- 10^-5

  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  conditions <- captureConditions(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(2020, 2020),
                                    harmonization = "absoluteChanges")
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  x <- result$x

  capTg <- 1200 * 10^-5 * 10^-3
  # the collapsed cell lands exactly on the cap instead of reaching 1e9
  expect_equal(as.vector(x["regone.1", 2025, "fertilizer"]), 1200)
  # its excess of almost 10 Tg is absorbed by its region partner
  expect_equal(as.vector(x["regone.2", 2025, "fertilizer"]),
               (10 + (10 - capTg)) / 100 * (10^9 / 10^6))
  # the partner in the other region is untouched
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer"]), 10)
  expect_true(max(x[, , "fertilizer"]) <= 1200)

  # tonnage and global totals stay intact and parallel to the (flat) input
  expect_equal(as.vector(dimSums(toolFertilizerTg(x[, 2025, "fertilizer"], landHarmonized), c(1, 2))), 21)
  expect_equal(as.vector(x["regone.1", 2020, "fertilizer"]), 1000)
})

test_that("absoluteChanges compensates lost fertilizer tonnage within the same region", {
  # regone.2 loses all its cropland in 2025, so its 10 Tg anchor is lost and
  # must be compensated by its own region: expected total 20 + 50 = 70 Tg over
  # an actual 60 Tg gives factor 7/6 -> 700 kg ha-1 yr-1 on the remaining cell,
  # while regtwo, which lost nothing, keeps its anchor rate of 100 untouched
  # (a global factor would inflate it to 114.29)
  cells <- c("regone.1", "regone.2", "regtwo.1")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput["regone.1", , "fertilizer.c3ann"] <- c(100, 100)
  xInput["regone.2", , "fertilizer.c3ann"] <- c(100, 150)
  xInput["regtwo.1", , "fertilizer.c3ann"] <- c(10, 10)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  # target rate of 100 kg ha-1 yr-1 on 100 Mha anchors 10 Tg per cell
  xTarget <- new.magpie(cells, years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = "c3ann", fill = 100)
  landHarmonized["regone.2", 2025, ] <- 0

  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  conditions <- captureConditions(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(2020, 2020),
                                    harmonization = "absoluteChanges")
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  x <- result$x

  # the lost 10 Tg are compensated by regone's own remaining cell, regtwo is
  # left at its anchor rate, the rates stay far below the 1200 cap
  expect_equal(as.vector(x["regone.1", 2025, "fertilizer"]), 700)
  expect_equal(as.vector(x["regone.2", 2025, "fertilizer"]), 0)
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer"]), 100)

  # each region's total equals its target total plus its own input delta
  regionTg <- dimSums(toolFertilizerTg(x[, 2025, "fertilizer"], landHarmonized), 1.2)
  expect_equal(as.vector(regionTg["regone", , ]), 20 + 50)
  expect_equal(as.vector(regionTg["regtwo", , ]), 10)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
})

test_that("absoluteChanges falls back to a global factor if a region lost all its cropland", {
  # regone expects its 10 Tg anchor in 2025 but lost all its cropland, so no
  # cell of the region can hold it: the compensation falls back to one global
  # factor, regone's warn check fires, the per-region check fails by design
  # while the global parallelism check stays green
  cells <- c("regone.1", "regtwo.1")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput["regone.1", , "fertilizer.c3ann"] <- c(100, 100)
  xInput["regtwo.1", , "fertilizer.c3ann"] <- c(10, 10)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = "c3ann", fill = 100)
  landHarmonized["regone.1", 2025, ] <- 0

  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  conditions <- captureConditions(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(2020, 2020),
                                    harmonization = "absoluteChanges")
  )
  # the warn check reports the unachievable region, the new per-region check
  # fails by design here, the global one stays green
  expect_true(any(grepl("Check failed.*cropland left to hold it", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional fertilizer totals are parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global fertilizer totals are parallel", conditions$conditions)))
  x <- result$x

  # the global factor 20/10 doubles regtwo's rate, regone holds nothing
  expect_equal(as.vector(x["regone.1", 2025, "fertilizer"]), 0)
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer"]), 200)
})

test_that("the global fallback is scoped to the offending crop group and year table", {
  # regtwo lost all c3ann cropland while still expecting its 10 Tg anchor, so
  # the c3ann 2025 table falls back to the global factor (30/10 -> regone.1 at
  # 300), but the c4ann table has no such unholdable region: regone lost its
  # regone.2 anchor to zero area and must be compensated by its own region
  # alone (20/10 -> 200), leaving regtwo's c4ann rate at its 100 anchor
  # (a global fallback for every table would give 150/150)
  cells <- c("regone.1", "regone.2", "regtwo.1")
  items <- c(nonlandItems, "fertilizer.c4ann")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = items, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput["regone.1", , c("fertilizer.c3ann", "fertilizer.c4ann")] <- 100
  xInput["regone.2", , c("fertilizer.c3ann", "fertilizer.c4ann")] <- 100
  xInput["regtwo.1", , c("fertilizer.c3ann", "fertilizer.c4ann")] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  # target rate of 100 kg ha-1 yr-1 on 100 Mha anchors 10 Tg per cell and crop
  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , c("fertilizer.c3ann", "fertilizer.c4ann")] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = c("c3ann", "c4ann"), fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = c("c3ann", "c4ann"), fill = 100)
  landHarmonized["regone.2", 2025, c("c3ann", "c4ann")] <- 0
  landHarmonized["regtwo.1", 2025, "c3ann"] <- 0

  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  conditions <- captureConditions(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(2020, 2020),
                                    harmonization = "absoluteChanges")
  )
  # the warn fires for regtwo's unholdable c3ann tonnes, the per-region check
  # fails only for the c3ann 2025 table, the global one stays green
  expect_true(any(grepl("Check failed.*cropland left to hold it", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional fertilizer totals are parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global fertilizer totals are parallel", conditions$conditions)))
  x <- result$x

  expect_equal(as.vector(x["regone.1", 2025, "fertilizer.c3ann"]), 300)
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer.c3ann"]), 0)
  expect_equal(as.vector(x["regone.1", 2025, "fertilizer.c4ann"]), 200)
  expect_equal(as.vector(x["regone.2", 2025, "fertilizer.c4ann"]), 0)
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer.c4ann"]), 100)

  # the c4ann regional totals stay exact despite the fallback in c3ann
  regionTg <- dimSums(toolFertilizerTg(x[, 2025, "fertilizer.c4ann"], landHarmonized), 1.2)
  expect_equal(as.vector(regionTg["regone", , ]), 20)
  expect_equal(as.vector(regionTg["regtwo", , ]), 10)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
})

test_that("stranded bioh of clusters with zero harvest area is redistributed", {
  # identical to the first test, except reg.one loses its harmonized harvest
  # area after 2020 while keeping positive bioh: rates x land cannot host it
  # where the land vanished, so the bioh of reg.one moves to
  # reg.two, the only cluster that still has harvest area
  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025, 2030), names = nonlandItems, fill = 0)
  xInput["reg.one", 2020, ] <- c(2000, 1500, 1, 1, 3, 1, 400)
  xInput["reg.one", 2025, ] <- c(500, 1600, 1, 1, 2, 2, 500)
  xInput["reg.one", 2030, ] <- c(500, 1600, 1, 1, 2, 2, 450)
  xInput["reg.two", 2020, ] <- c(1000, 1000, 1, 1, 1, 1, 500)
  xInput["reg.two", 2025, ] <- c(1200, 1000, 1, 1, 1, 1, 600)
  xInput["reg.two", 2030, ] <- c(1200, 1000, 1, 1, 1, 1, 550)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget["reg.one", 2010, ] <- c(1000, 1500, 1, 1, 2, 2, 500)
  xTarget["reg.one", 2020, ] <- c(1000, 1500, 1, 1, 2, 2, 500)
  xTarget["reg.two", 2010, ] <- c(1000, 1000, 1, 1, 1, 1, 400)
  xTarget["reg.two", 2020, ] <- c(1000, 1000, 1, 1, 1, 1, 400)
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020),
                           names = c("c3ann", "pastr"), fill = 100)
  landHarmonized <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020, 2025, 2030),
                               names = c("c3ann", "pastr"), fill = 100)
  landHarmonized["reg.two", 2025, "c3ann"] <- 0

  harvestArea <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020, 2025, 2030),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  harvestArea["reg.one", c(2025, 2030), ] <- 0
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetExtrapolated = xTarget,
             LandTargetExtrapolated = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  conditions <- captureConditions(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(2020, 2020),
                                    harmonization = "absoluteChanges")
  )
  # the total keeping check for bioh passes by construction, only the
  # implausible fertilizer rate of the unchanged fixture is reported
  expect_false(any(grepl("Check failed: Adapting bioh", conditions$conditions)))
  expect_true(any(grepl("Fertilizer application", conditions$conditions)))
  # the report covers the two stranded cluster years of reg.one (2025/2030)
  # with 1600 kg C yr-1 each, 3200 of the 16600 kg C yr-1 of all cluster years
  report <- grep("redistributed stranded bioh", conditions$conditions, value = TRUE)
  expect_equal(length(report), 1)
  expect_true(grepl("3200 kg C yr-1 \\(19\\.3%", report))
  expect_true(grepl("from 2 of 8 cluster years \\(25%\\)", report))
  x <- result$x

  # the stranded clusters end up exactly at zero
  expect_equal(as.vector(x["reg.one", c(2025, 2030), paste0("bioh.", c("primf", "secmf"))]),
               rep(0, 4))

  # reg.two hosts the harmonizer bioh of both clusters: reg.one keeps 1600 kg C
  # in 2025/2030 (primf clamped to 0, secmf 1600), which has no receiver in its
  # region and moves over the global pool to reg.two (1200 + 1000), normalization
  # then rescales the total of 3800 onto the calculated 1200/1000 split
  expect_equal(as.vector(x["reg.two", 2025, c("bioh.primf", "bioh.secmf")]),
               c(1200, 1000) * 3800 / 2200)
  expect_equal(as.vector(x["reg.two", 2030, c("bioh.primf", "bioh.secmf")]),
               c(1200, 1000) * 3800 / 2200)

  # global totals equal the pre-redistribution harmonizer totals
  expect_equal(as.numeric(dimSums(x[, 2025, "bioh"], c(1, 2, 3))), 1600 + 2200)
  expect_equal(as.numeric(dimSums(x[, 2030, "bioh"], c(1, 2, 3))), 1600 + 2200)

  # everything besides bioh and the harvest area stays as in the first test
  expect_equal(as.vector(x["reg.one", 2025, "fertilizer"]), 2900)
  expect_equal(as.vector(x["reg.two", 2025, "fertilizer"]), 0)
  expect_equal(as.vector(x["reg.one", 2030, "fertilizer"]), 1000)
  expect_equal(as.vector(x["reg.two", 2030, "fertilizer"]), 900)
  expect_equal(as.vector(x["reg.one", 2025, c("harvest_weight_type.roundwood",
                                              "harvest_weight_type.fuelwood")]), c(1, 3))
  expect_equal(as.vector(x["reg.one", 2025, "wood_harvest_area.primf"]), 0)
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
  expect_true(all(is.finite(x)))
  expect_true(all(x >= 0))
})

test_that("toolCapFertilizerKgPerHa redistributes within the region", {
  fert <- new.magpie(c("regone.1", "regone.2", "regtwo.1"), years = 2025,
                     names = "fertilizer.c3ann", fill = 0)
  fert["regone.1", , ] <- 1500
  fert["regone.2", , ] <- 100
  fert["regtwo.1", , ] <- 10
  getSets(fert) <- c("region", "id", "year", "category", "data")
  land <- new.magpie(c("regone.1", "regone.2", "regtwo.1"), years = 2025, names = "c3ann", fill = 100)

  capped <- toolCapFertilizerKgPerHa(fert, land)

  # regone's 30 Tg excess fits into the 110 Tg headroom of its own region
  # (cap 120 Tg per cell), so regtwo stays untouched
  expect_equal(as.vector(capped["regone.1", , ]), 1200)
  expect_equal(as.vector(capped["regone.2", , ]), 400)
  expect_equal(as.vector(capped["regtwo.1", , ]), 10)
  expect_equal(as.vector(dimSums(toolFertilizerTg(capped, land), c(1, 2))), 150 + 10 + 1)
})

test_that("toolCapFertilizerKgPerHa falls back to the global pool and skips infeasible cases", {
  fert <- new.magpie(c("regone.1", "regtwo.1"), years = 2025, names = "fertilizer.c3ann", fill = 0)
  fert["regone.1", , ] <- 1500
  fert["regtwo.1", , ] <- 100
  getSets(fert) <- c("region", "id", "year", "category", "data")
  land <- new.magpie(c("regone.1", "regtwo.1"), years = 2025, names = "c3ann", fill = 100)

  # regone has no headroom of its own, the pass runs over the global pool:
  # its 30 Tg excess fits into regtwo's 110 Tg spare
  capped <- toolCapFertilizerKgPerHa(fert, land)
  expect_equal(as.vector(capped["regone.1", , ]), 1200)
  expect_equal(as.vector(capped["regtwo.1", , ]), 400)

  # if not even the global pool has enough headroom, tonnage loss would be
  # unavoidable, so the input is returned unchanged for toolCheckFertilizer to
  # report
  land["regtwo.1", , ] <- 10
  expect_equal(as.vector(toolCapFertilizerKgPerHa(fert, land)), as.vector(fert))
})
