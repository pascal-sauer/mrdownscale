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

  # the parallel target 4500 + (4300 - 5500) = 3300 of the merged region
  # (dimSums reads the dots in reg.one/reg.two as cell separators, merging both
  # cells into one region) is allocated directly onto biohCalculated: the
  # clamped kg C per Mha rate of 0 for reg.one primf is floored to the smallest
  # positive rate (1000), giving weights 1000/1600 in reg.one and 1200/1000 in
  # reg.two, so every cluster holds its weight share of the 4800 of the region
  expect_equal(as.vector(x["reg.one", 2025, c("bioh.primf", "bioh.secmf")]),
               3300 * c(1000, 1600) / 4800)
  expect_equal(as.vector(x["reg.two", 2025, c("bioh.primf", "bioh.secmf")]),
               3300 * c(1200, 1000) / 4800)
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
  # regone targets its 10 Tg anchor in 2025 but lost all its cropland, so no
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

test_that("zeroing cropland-less fertilizer tonnages matches clusters by name", {
  # the mocked nonland input lists its clusters in the reversed order of the
  # nonland target: the projected tonnages follow the target's item order, the
  # cropland mask the input's, so a positional mask zeroes the wrong cluster
  xInput <- new.magpie(c("regtwo.1", "regone.1"), years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput["regone.1", , "fertilizer.c3ann"] <- c(10, 10)
  xInput["regtwo.1", , "fertilizer.c3ann"] <- c(10, 10)
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  # target rate of 100 kg ha-1 yr-1 on 100 Mha anchors 10 Tg per cluster
  xTarget <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020, 2025),
                               names = "c3ann", fill = 100)
  # regtwo.1 loses all its cropland in 2025, so only its 10 Tg anchor must be zeroed
  landHarmonized["regtwo.1", 2025, ] <- 0

  harvestArea <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020, 2025),
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
  # regtwo lost the land holding its tonnes, its warn fires, while the global
  # parallelism check must stay green
  expect_true(any(grepl("Check failed.*cropland left to hold it", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global fertilizer totals are parallel", conditions$conditions)))
  x <- result$x

  # only regtwo's tonnes are zeroed, its unholdable target falls back globally
  # onto regone, doubling its 10 Tg anchor to a 200 kg ha-1 yr-1 rate
  expect_equal(as.vector(x["regone.1", 2025, "fertilizer"]), 200)
  expect_equal(as.vector(x["regtwo.1", 2025, "fertilizer"]), 0)
})

test_that("clusters with zero harvest area host no bioh for absoluteChanges", {
  # identical to the first test, except reg.one loses its harmonized harvest
  # area after 2020 while keeping positive bioh: under absoluteChanges the
  # regional cover target is allocated onto biohCalculated, so the
  # weightless reg.one clusters host nothing and the whole target of the
  # merged region lands on reg.two, the only cluster with harvest area
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
  # no redistribution happens for absoluteChanges, the allocation gives the
  # weightless clusters zero by construction, no bioh check fails and only the
  # implausible fertilizer rate of the unchanged fixture is reported
  expect_false(any(grepl("Check failed.*bioh", conditions$conditions)))
  expect_true(any(grepl("Fertilizer application", conditions$conditions)))
  x <- result$x

  # the weightless clusters end up exactly at zero
  expect_equal(as.vector(x["reg.one", c(2025, 2030), paste0("bioh.", c("primf", "secmf"))]),
               rep(0, 4))

  # reg.two hosts the whole parallel target of the merged region (dots
  # read as cell separators merge the cells into one region totalling
  # 4500 + (4300 - 5500) = 3300) on the only weight left, its own
  # recalculated 1200/1000 split, so 3300 * c(1200, 1000) / 2200
  expect_equal(as.vector(x["reg.two", 2025, c("bioh.primf", "bioh.secmf")]),
               c(1200, 1000) * 3300 / 2200)
  expect_equal(as.vector(x["reg.two", 2030, c("bioh.primf", "bioh.secmf")]),
               c(1200, 1000) * 3300 / 2200)

  # global totals equal the regional rescaled parallel total 4500 - 1200
  expect_equal(as.numeric(dimSums(x[, 2025, "bioh"], c(1, 2, 3))), 4500 - 1200)
  expect_equal(as.numeric(dimSums(x[, 2030, "bioh"], c(1, 2, 3))), 4500 - 1200)

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

test_that("absoluteChanges keeps regional and global bioh totals parallel to growing input", {
  # regone.1's input bioh grows while land and harvest area stay flat: nothing
  # is clamped, so every region's bioh total equals its target total of the
  # harmonization year plus its own input delta, also globally
  cells <- c("regone.1", "regone.2", "regtwo.1")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.1", 2025, c("bioh.primf", "bioh.secmf")] <- c(1500, 1200)
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
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
  expect_false(any(grepl("Check failed.*bioh", conditions$conditions)))
  x <- result$x

  # regone: target 4000 + delta 700, regtwo: target 2000 + delta 0
  biohRegion <- dimSums(x[, 2025, "bioh"], c(1.2, 3))
  expect_equal(as.vector(biohRegion["regone", , ]), 4000 + 700)
  expect_equal(as.vector(biohRegion["regtwo", , ]), 2000)
  expect_equal(as.numeric(dimSums(x[, 2025, "bioh"], c(1, 2, 3))), 6000 + 700)

  # with parallelism already holding per cell, the cells keep their values
  expect_equal(as.vector(x["regone.1", 2025, c("bioh.primf", "bioh.secmf")]), c(1500, 1200))
  expect_equal(as.vector(x["regone.2", 2025, c("bioh.primf", "bioh.secmf")]), c(1000, 1000))

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
})

test_that("absoluteChanges compensates clamped bioh within the region", {
  # regone.1's input bioh collapses far below the target, clamping its primf to
  # zero: the region's parallel target 4000 - 2500 + 500 = 2000 is
  # allocated onto the recalculated bioh (primf floored to the smallest positive
  # rate 500), so both cells share it proportional to their recalculated
  # 1000 and 2500, the other cell absorbing the difference
  cells <- c("regone.1", "regone.2", "regtwo.1")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.1", 2020, c("bioh.primf", "bioh.secmf")] <- c(3000, 1000)
  xInput["regone.1", 2025, c("bioh.primf", "bioh.secmf")] <- c(1000, 500)
  xInput["regone.2", 2025, "bioh.primf"] <- 1500
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
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
  expect_false(any(grepl("Check failed.*bioh", conditions$conditions)))
  x <- result$x

  biohRegion <- dimSums(x[, 2025, "bioh"], c(1.2, 3))
  expect_equal(as.vector(biohRegion["regone", , ]), 4000 - 2500 + 500)
  expect_equal(as.vector(biohRegion["regtwo", , ]), 2000)
  expect_equal(as.numeric(dimSums(x[, 2025, "bioh"], c(1, 2, 3))), 2000 + 2000)

  # both cells of regone share the target 2000 proportional to their
  # recalculated weights, regone.1 on its floored 500/500, regone.2 on 1500/1000
  expect_equal(as.vector(x["regone.1", 2025, c("bioh.primf", "bioh.secmf")]),
               2000 * c(500, 500) / 3500)
  expect_equal(as.vector(x["regone.2", 2025, c("bioh.primf", "bioh.secmf")]),
               2000 * c(1500, 1000) / 3500)
  expect_equal(as.vector(x["regtwo.1", 2025, c("bioh.primf", "bioh.secmf")]), c(1000, 1000))

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
})

test_that("absoluteChanges allocates regional bioh proportional to recalculated bioh", {
  # within a region-cover the bioh of the clusters is exactly proportional to
  # biohCalculated (harmonized kg C per Mha times harmonized harvest area):
  # regone.1 recalculates to 2 x 7200, regone.2 to 2 x 1200 and regone.3 to
  # zero, so of the region target 9200 regone.1 holds 6/7, regone.2 1/7
  # and the capacity-less regone.3 exactly nothing, the ratio of the hosting
  # clusters being their recalculated ratio 6:1
  cells <- c("regone.1", "regone.2", "regone.3")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput["regone.1", , c("bioh.primf", "bioh.secmf")] <- 2000
  xInput["regone.1", 2025, c("bioh.primf", "bioh.secmf")] <- 2400
  xInput["regone.2", , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.2", 2025, c("bioh.primf", "bioh.secmf")] <- 1200
  xInput["regone.3", , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget["regone.1", , c("bioh.primf", "bioh.secmf")] <- 2000
  xTarget["regone.2", , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget["regone.3", , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = "c3ann", fill = 100)

  # the harmonized harvest area of regone.1 is tripled, regone.3 loses all of
  # its harvest area, so the recalculated bioh weights are 14400, 2400 and 0
  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  harvestArea["regone.1", 2025, ] <- 3
  harvestArea["regone.3", 2025, ] <- 0
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

  # the region target 8000 + (9200 - 8000) = 9200 is allocated onto the
  # recalculated 14400 : 2400 : 0, so regone.1 and regone.2 hold 6/7 and 1/7
  # of it, split onto their items by their recalculated item weights
  expect_equal(as.vector(x["regone.1", 2025, c("bioh.primf", "bioh.secmf")]),
               rep(9200 * 7200 / 16800, 2))
  expect_equal(as.vector(x["regone.2", 2025, c("bioh.primf", "bioh.secmf")]),
               rep(9200 * 1200 / 16800, 2))
  expect_equal(as.vector(x["regone.3", 2025, c("bioh.primf", "bioh.secmf")]), c(0, 0))
  expect_equal(as.numeric(sum(x["regone.1", 2025, c("bioh.primf", "bioh.secmf")]) /
                            sum(x["regone.2", 2025, c("bioh.primf", "bioh.secmf")])), 6)
  expect_equal(as.numeric(dimSums(x[, 2025, "bioh"], c(1, 2, 3))), 9200)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)
})

test_that("absoluteChanges leaves stranded bioh unallocated when no global room is left, failing both checks", {
  # regone.1 keeps positive input bioh but loses all its harvest area after 2020
  # and has no receiver in its region: its parallel target of 2500 - 1400
  # = 1100 kg C yr-1 finds no global room, because regtwo.1's own share of
  # 2000 already fills its recalculated weight, so the stranded 1100 stays
  # unallocated and the regional and global checks fail by design
  xInput <- new.magpie(c("regone.1", "regtwo.1"), years = c(2020, 2025, 2030), names = nonlandItems, fill = 0)
  xInput["regone.1", 2020, c("bioh.primf", "bioh.secmf")] <- c(2000, 1500)
  xInput["regone.1", c(2025, 2030), "bioh.primf"] <- 500
  xInput["regone.1", c(2025, 2030), "bioh.secmf"] <- 1600
  xInput["regtwo.1", , "bioh.primf"] <- 1000
  xInput["regtwo.1", , "bioh.secmf"] <- 1000
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020), names = nonlandItems, fill = 0)
  xTarget["regone.1", , "bioh.primf"] <- 1000
  xTarget["regone.1", , "bioh.secmf"] <- 1500
  xTarget["regtwo.1", , "bioh.primf"] <- 1000
  xTarget["regtwo.1", , "bioh.secmf"] <- 1000
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020, 2025, 2030),
                               names = "c3ann", fill = 100)

  harvestArea <- new.magpie(c("regone.1", "regtwo.1"), years = c(2010, 2020, 2025, 2030),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  harvestArea["regone.1", c(2025, 2030), ] <- 0
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
  expect_true(any(grepl("Check failed.*Regional forest bioh totals are parallel", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Global bioh totals are parallel", conditions$conditions)))
  x <- result$x

  # regone.1 hosts nothing, regtwo.1 hosts only its own 2000, the stranded
  # 1100 finds no remaining room and is lost
  biohRegion <- dimSums(x[, c(2025, 2030), "bioh"], c(1.2, 3))
  expect_equal(as.numeric(biohRegion["regone", , ]), rep(0, 2))
  expect_equal(as.numeric(biohRegion["regtwo", , ]), rep(2000, 2))
  expect_equal(as.numeric(dimSums(x[, 2025, "bioh"], c(1, 2, 3))), 2000)
  expect_equal(as.numeric(dimSums(x[, 2030, "bioh"], c(1, 2, 3))), 2000)
})

test_that("absoluteChanges keeps regional forest and other land cover totals of bioh parallel to the input", {
  # bioh item covers: forest = primf, secmf, other = secnf. regone.1 grows its
  # forest bioh while its other land bioh collapses from 500 to 1 far below its
  # target anchor of 300, so the raw absolute change clamps to zero inside one
  # cover. The per-cell renormalization onto kg C per Mha times harvest area
  # then re-splits the cell totals between the covers, which the final
  # cover-granular rescale must undo, keeping each region's cover total at its
  # floored target total plus regional input cover delta
  cells <- c("regone.1", "regone.2", "regtwo.1")
  coverItems <- c("bioh.primf", "bioh.secmf", "bioh.secnf")
  items <- c(nonlandItems, "bioh.secnf", "wood_harvest_area.secnf")

  xInput <- new.magpie(cells, years = c(2020, 2025), names = items, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.1", , "bioh.secnf"] <- 500
  xInput["regone.1", 2025, "bioh.secnf"] <- 1
  xInput["regone.1", 2025, "bioh.secmf"] <- 3000
  xInput["regone.2", , "bioh.secnf"] <- 100
  xInput["regone.2", 2025, "bioh.secnf"] <- 300
  xInput["regtwo.1", , "bioh.secnf"] <- 100
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
               "wood_harvest_area.secnf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget["regone.1", , "bioh.secnf"] <- 300
  xTarget["regone.2", , "bioh.secnf"] <- 100
  xTarget["regtwo.1", , "bioh.secnf"] <- 100
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                "wood_harvest_area.secnf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025), names = "c3ann", fill = 100)

  # both covers keep 1 Mha yr-1 of harvest area everywhere
  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                                      "wood_harvest_area.secnf"), fill = 1)
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

  # floored target cover totals plus regional input cover deltas: forest 4000 +
  # 2000 = 6000 in regone and 2000 in regtwo, other land 400 - 299 = 101 in
  # regone and 100 in regtwo, also globally
  forest <- dimSums(x[, 2025, c("bioh.primf", "bioh.secmf")], c(1.2, 3))
  other <- dimSums(x[, 2025, "bioh.secnf"], c(1.2, 3))
  expect_equal(as.vector(forest["regone", , ]), 6000, tolerance = 10^-6)
  expect_equal(as.vector(forest["regtwo", , ]), 2000, tolerance = 10^-6)
  expect_equal(as.vector(other["regone", , ]), 101, tolerance = 10^-6)
  expect_equal(as.vector(other["regtwo", , ]), 100, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(forest, 1)[, 2025, ]), 8000, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(other, 1)[, 2025, ]), 201, tolerance = 10^-6)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), coverItems] - xTarget[, c(2010, 2020), coverItems])) < 10^-5)
})

test_that("absoluteChanges allocates the global other land cover target first when a regional one is negative", {
  # regone's input other land bioh collapses from 300 to 1 while its target
  # anchor in the harmonization year is only 100, so its un-floored other land
  # target turns negative while regtwo grows: flooring the regional
  # targets per region would drop regone's negative delta and lift the
  # global other land total above the raw global target, so the global
  # first allocation gives regone nothing and scales regtwo down to keep the
  # global cover total exactly parallel to the un-floored target
  cells <- c("regone.1", "regtwo.1")
  coverItems <- c("bioh.primf", "bioh.secmf", "bioh.secnf")
  items <- c(nonlandItems, "bioh.secnf", "wood_harvest_area.secnf")

  # the forest bioh stays flat at its target anchor in both regions
  xInput <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.1", , "bioh.secnf"] <- c(300, 1, 1)
  xInput["regtwo.1", , "bioh.secnf"] <- c(100, 400, 600)
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
               "wood_harvest_area.secnf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget["regone.1", , "bioh.secnf"] <- 100
  xTarget["regtwo.1", , "bioh.secnf"] <- 100
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                "wood_harvest_area.secnf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025, 2030), names = "c3ann", fill = 100)

  # both covers keep 1 Mha yr-1 of harvest area everywhere
  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025, 2030),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                                      "wood_harvest_area.secnf"), fill = 1)
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

  # the un-floored other land targets are regone 100 - 299 = -199 in both
  # future years and regtwo 400 in 2025 and 600 in 2030, giving the raw global
  # targets 201 and 401
  other <- dimSums(x[, c(2025, 2030), "bioh.secnf"], c(1.2, 3))

  # the global other land total equals the un-floored raw global target
  expect_equal(as.vector(dimSums(other, 1)[, 2025, ]), 201, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(other, 1)[, 2030, ]), 401, tolerance = 10^-9)

  # the loss-making region gets nothing of the cover
  expect_equal(as.vector(other["regone", 2025, ]), 0)
  expect_equal(as.vector(other["regone", 2030, ]), 0)

  # the growing region carries the whole global target, scaled below its own
  expect_equal(as.vector(other["regtwo", 2025, ]), 201, tolerance = 10^-9)
  expect_equal(as.vector(other["regtwo", 2030, ]), 401, tolerance = 10^-9)

  # the forest cover with its unharmed targets stays exactly the same
  forest <- dimSums(x[, c(2025, 2030), c("bioh.primf", "bioh.secmf")], c(1.2, 3))
  expect_equal(as.vector(forest["regone", 2025, ]), 2000, tolerance = 10^-9)
  expect_equal(as.vector(forest["regtwo", 2030, ]), 2000, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forest, 1)[, 2025, ]), 4000, tolerance = 10^-9)
  expect_equal(as.vector(dimSums(forest, 1)[, 2030, ]), 4000, tolerance = 10^-9)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), coverItems] - xTarget[, c(2010, 2020), coverItems])) < 10^-5)
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges charges a negative global other land cover target to the forest cover", {
  # both regions' input other land bioh collapses far below their target anchors
  # so the global other land target turns negative while the flat global
  # forest target outweighs the loss: the floored other land total must be
  # charged to the forest cover so the global total over all covers stays
  # parallel to the input and the positive forest regions are uniformly scaled
  # by the balancing factor 3702 / 4000
  cells <- c("regone.1", "regtwo.1")
  coverItems <- c("bioh.primf", "bioh.secmf", "bioh.secnf")
  items <- c(nonlandItems, "bioh.secnf", "wood_harvest_area.secnf")

  # the forest bioh stays flat at its target anchor in both regions
  xInput <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.1", , "bioh.secnf"] <- c(300, 1, 1)
  xInput["regtwo.1", , "bioh.secnf"] <- c(200, 1, 1)
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
               "wood_harvest_area.secnf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget["regone.1", , "bioh.secnf"] <- 100
  xTarget["regtwo.1", , "bioh.secnf"] <- 100
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                "wood_harvest_area.secnf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025, 2030), names = "c3ann", fill = 100)

  # both covers keep 1 Mha yr-1 of harvest area everywhere
  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025, 2030),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                                      "wood_harvest_area.secnf"), fill = 1)
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

  # the un-floored other land targets are regone 100 - 299 = -199 and
  # regtwo 100 - 199 = -99, giving the global targets -298 for other land
  # and 4000 for forest in both future years
  other <- dimSums(x[, c(2025, 2030), "bioh.secnf"], c(1.2, 3))
  forest <- dimSums(x[, c(2025, 2030), c("bioh.primf", "bioh.secmf")], c(1.2, 3))

  # the floored other land cover gives up its whole global total
  expect_equal(as.vector(dimSums(other, 1)[, 2025, ]), 0)
  expect_equal(as.vector(dimSums(other, 1)[, 2030, ]), 0)

  # the forest cover absorbs the charged loss: its global total is the summed
  # un-floored global target 4000 - 298, not its own global target
  expect_equal(as.vector(dimSums(forest, 1)[, 2025, ]), 3702, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(forest, 1)[, 2030, ]), 3702, tolerance = 10^-6)

  # the global total over both covers equals the un-floored target 3702
  expect_equal(as.vector(dimSums(forest + other, 1)[, 2025, ]), 3702, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(forest + other, 1)[, 2030, ]), 3702, tolerance = 10^-6)

  # both positive forest regions are scaled by the balancing factor 3702 / 4000
  expect_equal(as.vector(forest["regone", 2025, ]), 2000 * (3702 / 4000), tolerance = 10^-6)
  expect_equal(as.vector(forest["regtwo", 2025, ]), 2000 * (3702 / 4000), tolerance = 10^-6)
  expect_equal(as.vector(forest["regone", 2030, ]), 2000 * (3702 / 4000), tolerance = 10^-6)
  expect_equal(as.vector(forest["regtwo", 2030, ]), 2000 * (3702 / 4000), tolerance = 10^-6)

  # before and in the harmonization year the target is still reproduced exactly
  expect_true(max(abs(x[, c(2010, 2020), coverItems] - xTarget[, c(2010, 2020), coverItems])) < 10^-5)
  expect_true(min(x) >= 0)
})

test_that("absoluteChanges cross-cover balancing is a no-op without a negative global cover target", {
  # regone's other land target collapses below zero but regtwo's growth
  # keeps the global other land target positive: the balancing factor is
  # exactly 1 and both cover global totals stay parallel as before
  cells <- c("regone.1", "regtwo.1")
  coverItems <- c("bioh.primf", "bioh.secmf", "bioh.secnf")
  items <- c(nonlandItems, "bioh.secnf", "wood_harvest_area.secnf")

  xInput <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xInput["regone.1", , "bioh.secnf"] <- c(300, 1, 1)
  xInput["regtwo.1", , "bioh.secnf"] <- c(100, 500, 500)
  xInput[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
               "wood_harvest_area.secnf")] <- 1
  xInput[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xInput[, , "fertilizer.c3ann"] <- 100
  attr(xInput, "geometry") <- "point"
  attr(xInput, "crs") <- "EPSG:4326"
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , c("bioh.primf", "bioh.secmf")] <- 1000
  xTarget["regone.1", , "bioh.secnf"] <- 100
  xTarget["regtwo.1", , "bioh.secnf"] <- 100
  xTarget[, , c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                "wood_harvest_area.secnf")] <- 1
  xTarget[, , c("harvest_weight_type.roundwood", "harvest_weight_type.fuelwood")] <- 1
  xTarget[, , "fertilizer.c3ann"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  landTarget <- new.magpie(cells, years = c(2010, 2020), names = "c3ann", fill = 100)
  landHarmonized <- new.magpie(cells, years = c(2010, 2020, 2025, 2030), names = "c3ann", fill = 100)

  harvestArea <- new.magpie(cells, years = c(2010, 2020, 2025, 2030),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf",
                                      "wood_harvest_area.secnf"), fill = 1)
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

  # the un-floored other land targets are regone -199 and regtwo 500, the
  # global targets 4000 for forest and 301 for other land, both positive
  other <- dimSums(x[, c(2025, 2030), "bioh.secnf"], c(1.2, 3))
  forest <- dimSums(x[, c(2025, 2030), c("bioh.primf", "bioh.secmf")], c(1.2, 3))

  # both cover global totals equal their own un-floored global targets
  expect_equal(as.vector(dimSums(forest, 1)[, 2025, ]), 4000, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(forest, 1)[, 2030, ]), 4000, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(other, 1)[, 2025, ]), 301, tolerance = 10^-6)
  expect_equal(as.vector(dimSums(other, 1)[, 2030, ]), 301, tolerance = 10^-6)

  # the loss-making region stays at zero, the growing region keeps its floored
  # and scaled total, both forest regions stay parallel
  expect_equal(as.vector(other["regone", 2025, ]), 0)
  expect_equal(as.vector(other["regtwo", 2025, ]), 301, tolerance = 10^-6)
  expect_equal(as.vector(forest["regone", 2025, ]), 2000, tolerance = 10^-6)
  expect_equal(as.vector(forest["regtwo", 2025, ]), 2000, tolerance = 10^-6)

  expect_true(max(abs(x[, c(2010, 2020), coverItems] - xTarget[, c(2010, 2020), coverItems])) < 10^-5)
  expect_true(min(x) >= 0)
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
  # its 30 Tg excess fits into regtwo's 110 Tg spare, reported as a warning
  conditions <- captureConditions(capped <- toolCapFertilizerKgPerHa(fert, land))
  expect_true(any(grepl("Check failed.*Regions absorb their own fertilizer", conditions$conditions)))
  expect_equal(as.vector(capped["regone.1", , ]), 1200)
  expect_equal(as.vector(capped["regtwo.1", , ]), 400)

  # if not even the global pool has enough headroom, tonnage loss would be
  # unavoidable, so the input is returned unchanged for toolCheckFertilizer to
  # report
  land["regtwo.1", , ] <- 10
  expect_equal(as.vector(toolCapFertilizerKgPerHa(fert, land)), as.vector(fert))
})
