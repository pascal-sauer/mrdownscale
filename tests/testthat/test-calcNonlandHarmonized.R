nonlandItems <- c("bioh.primf", "bioh.secmf", "wood_harvest_area.primf", "wood_harvest_area.secmf",
                  "harvest_weight_type.roundwood", "harvest_weight_type.fuelwood", "fertilizer.c3ann")

test_that("calcNonlandHarmonized works with harmonization = absoluteChanges", {
  harmonizationYear <- 2020

  # real MAgPIE/LUH3 source data is not available in CI, so the source data
  # functions are mocked with synthetic data
  # for reg.one the kg C per Mha of primf drops in the input from 2000 to 500,
  # the raw absolute change is negative, so the harmonizer clamps it to zero
  # and it is later floored to the smallest positive rate
  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025), names = nonlandItems, fill = 0)
  xInput["reg.one", 2020, ] <- c(2000, 1500, 1, 1, 3, 1, 400)
  xInput["reg.one", 2025, ] <- c(500, 1600, 1, 1, 2, 2, 500)
  xInput["reg.two", 2020, ] <- c(1000, 1000, 1, 1, 1, 1, 500)
  xInput["reg.two", 2025, ] <- c(1200, 1000, 1, 1, 1, 1, 600)
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
  landHarmonized <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020, 2025),
                               names = c("c3ann", "pastr"), fill = 100)
  # reg.one raw fertilizer of 110 Tg divided by 100 Mha gives 1100 kg ha-1 yr-1,
  # reg.two has no cropland left in 2025, turning its raw fertilizer of 140 Tg
  # into an infinite rate, which must be set to zero
  landHarmonized["reg.two", 2025, "c3ann"] <- 0

  # wood harvest area stays at 1 Mha yr-1 everywhere
  harvestArea <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020, 2025),
                            names = c("wood_harvest_area.primf", "wood_harvest_area.secmf"), fill = 1)
  getSets(harvestArea) <- c("region", "id", "year", "category", "data")

  # temporarily replaces calcOutput in the package namespace, so the internal
  # source data calls return the synthetic data above and anything else errors
  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetLowRes = xTarget,
             LandTargetLowRes = landTarget,
             LandHarmonized = landHarmonized,
             WoodHarvestAreaHarmonized = harvestArea,
             # absoluteChanges only uses the target data up to the harmonization
             # year, so extrapolated targets must never be requested
             NonlandTargetExtrapolated = stop("\"", name, "\" must not be used for absoluteChanges"),
             LandTargetExtrapolated = stop("\"", name, "\" must not be used for absoluteChanges"),
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  expect_no_warning(
    result <- calcNonlandHarmonized(input = "magpie", target = "luh3",
                                    harmonizationPeriod = c(harmonizationYear, harmonizationYear),
                                    harmonization = "absoluteChanges")
  )
  x <- result$x

  expect_false(is.null(attr(x, "geometry")))
  expect_false(is.null(attr(x, "crs")))
  expect_true(setequal(getItems(x, 3), nonlandItems))
  expect_equal(getYears(x, as.integer = TRUE), c(2010, 2020, 2025))

  # before and in the harmonization year target data is returned unchanged
  expect_true(max(abs(x[, c(2010, 2020), ] - xTarget[, c(2010, 2020), ])) < 10^-5)

  # infinite fertilizer rates are replaced with zero, reg.one reproduces the
  # absolute changes of the input fertilizer applied to the target fertilizer
  # of the harmonization year: (50 + (500 - 400)) Tg / 100 Mha = 1500 kg ha-1 yr-1
  expect_equal(as.vector(x["reg.one", 2025, "fertilizer"]), 1500)
  expect_equal(as.vector(x["reg.two", 2025, "fertilizer"]), 0)
  expect_true(all(is.finite(x)))

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
             NonlandTargetLowRes = xTarget,
             LandTargetLowRes = land[, c(2010, 2020), ],
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
