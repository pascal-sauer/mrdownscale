items <- c("primf", "primn", "secdf", "secdn", "urban", "pastr", "c3ann")

test_that("calcLandTargetExtrapolated returns target unchanged if harmonization start and end year are equal", {
  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025), names = items, fill = 0)
  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = items, fill = 10)

  # if the guard does not return early, extrapolation, normalization or the
  # wood harvest area calculation would run and hit one of the stop() branches
  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandInputRecategorized = xInput,
             LandTargetLowRes = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    }
  )

  result <- calcLandTargetExtrapolated(input = "magpie", target = "luh3",
                                       harmonizationPeriod = c(2020, 2020))

  expect_identical(result$x, xTarget)
  expect_null(result$woodHarvestArea)
  expect_equal(result$unit, "Mha")
  expect_false(result$isocountries)
})

test_that("calcLandTargetExtrapolatedCore returns target unchanged if harmonization start and end year are equal", {
  xInput <- new.magpie(c("reg.one", "reg.two"), years = c(2020, 2025), names = items, fill = 0)
  xTarget <- new.magpie(c("reg.one", "reg.two"), years = c(2010, 2020), names = items, fill = 10)

  local_mocked_bindings(
    calcOutput = function(name, ...) {
      switch(name,
             LandInputRecategorized = xInput,
             LandTargetLowRes = xTarget,
             stop("unexpected calcOutput call for \"", name, "\""))
    },
    toolExtrapolate = function(x, years) stop("toolExtrapolate must not be used for equal harmonization years")
  )

  result <- calcLandTargetExtrapolatedCore(input = "magpie", target = "luh3",
                                           harmonizationPeriod = c(2020, 2020))

  expect_identical(result$x, xTarget)
})
