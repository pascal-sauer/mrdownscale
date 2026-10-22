nonlandItems <- c("bioh.primf", "bioh.secmf", "wood_harvest_area.primf", "wood_harvest_area.secmf",
                  "harvest_weight_type.roundwood", "harvest_weight_type.fuelwood", "fertilizer.c3ann")

test_that("calcNonlandTargetExtrapolated returns target unchanged if harmonization start and end year are equal", {
  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025), names = nonlandItems, fill = 1)
  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = nonlandItems, fill = 1)

  # if the guard does not return early, extrapolation would run and hit one of
  # the stop() branches (toolExtrapolate also fails on an empty year vector)
  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             NonlandInputRecategorized = xInput,
             NonlandTargetLowRes = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  result <- calcNonlandTargetExtrapolated(input = "magpie", target = "luh3",
                                          harmonizationPeriod = c(2020, 2020))

  expect_identical(result$x, xTarget)
  expect_true(grepl("fertilizer: kg ha-1 yr-1", result$unit))
  expect_false(result$isocountries)
})
