items <- c("primf", "secdf", "pltns", "primn", "secdn", "urc")

makeMagpie <- function(years, values) {
  # values: named list of numeric vectors (one per item, same length as years)
  arr <- array(as.numeric(unlist(lapply(items, function(i) values[[i]]))),
               dim = c(1, length(years), length(items)),
               dimnames = list("AFR", paste0("y", years), items))
  x <- as.magpie(arr, spatial = 1, temporal = 2)
  getSets(x) <- c("region", "year", "data")
  return(x)
}

test_that("core formula: target at harmonization year plus absolute changes of input", {
  target <- makeMagpie(c(2000, 2010, 2020), list(
    primf = c(105, 102, 100), secdf = c(48, 49, 50), pltns = c(5, 5, 5),
    primn = c(82, 81, 80), secdn = c(35, 38, 40), urc = c(30, 30, 30)))
  input <- makeMagpie(c(2020, 2025, 2030), list(
    primf = c(100, 100, 100), secdf = c(50, 52, 54), pltns = c(5, 5, 5),
    primn = c(80, 80, 80), secdn = c(40, 40, 40), urc = c(30, 28, 26)))

  expect_silent({
    out <- toolHarmonizeAbsoluteChanges(input, target, harmonizationPeriod = 2020)
  })

  expect_identical(getYears(out, as.integer = TRUE), c(2000L, 2010L, 2020L, 2025L, 2030L))
  expect_equal(as.vector(out[, 2020, ]), as.vector(target[, 2020, ]))
  expect_equal(
    as.vector(out[, 2025, ]),
    as.vector(target[, 2020, ] + (input[, 2025, ] - setYears(input[, 2020, ], NULL))))
  # total area conserved for every year
  expect_equal(as.vector(dimSums(out, 3)), rep(305, 5), tolerance = 10^-9)
  # passthrough of target before harmonization year
  expect_equal(as.vector(out[, c(2000, 2010), ]), as.vector(target[, c(2000, 2010), ]))
})

test_that("negative forest is funded from primn/secdn proportional to their shares", {
  target <- makeMagpie(c(2010, 2020), list(
    primf = c(100, 100), secdf = c(50, 50), pltns = c(0, 0),
    primn = c(60, 60), secdn = c(40, 40), urc = c(30, 30)))
  # input starts 20 Mha above target for secdf, drops to 0 -> out would be -20
  input <- makeMagpie(c(2020, 2025), list(
    primf = c(100, 100), secdf = c(70, 0), pltns = c(0, 0),
    primn = c(60, 60), secdn = c(40, 40), urc = c(10, 80)))

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(input, target, harmonizationPeriod = 2020)
  }, "1 negative forest cells \\(5\\.6 % of all cells\\).*min = -20.*mean = -20.*median = -20")

  expect_true(all(out >= 0))
  # shortfall of 20 deducted from primn (12) and secdn (8) proportional to their shares
  expect_equal(as.vector(out[, 2025, c("primn", "secdn")]), c(48, 32))
  expect_equal(as.vector(dimSums(out[, 2025, c("primn", "secdn")], 3)),
               as.vector(dimSums(target[, 2020, c("primn", "secdn")], 3)) - 20)
  # total conserved
  expect_equal(as.vector(dimSums(out, 3)), rep(280, 3), tolerance = 10^-9)
  # primf untouched by non-prim scaling
  expect_equal(as.vector(out[, 2025, "primf"]), 100)
})

test_that("other land insufficient -> clamp plus non-prim scaling", {
  target <- makeMagpie(c(2010, 2020), list(
    primf = c(100, 100), secdf = c(50, 50), pltns = c(0, 0),
    primn = c(10, 10), secdn = c(10, 10), urc = c(110, 110)))
  # shortfall of 40 vs only 20 Mha available in primn + secdn
  input <- makeMagpie(c(2020, 2025), list(
    primf = c(100, 100), secdf = c(90, 0), pltns = c(0, 0),
    primn = c(10, 10), secdn = c(10, 10), urc = c(90, 180)))

  expect_warning({
    out <- toolHarmonizeAbsoluteChanges(input, target, harmonizationPeriod = 2020)
  }, "negative forest cells")

  expect_true(all(out >= 0))
  # total equals target total at harmonization year for every year
  expect_equal(as.vector(dimSums(out, 3)), rep(280, 3), tolerance = 10^-9)
  # primf/primn untouched by non-prim scaling (primn fully spent funding secdf)
  expect_equal(as.vector(out[, 2025, "primf"]), 100)
  expect_equal(as.vector(out[, 2025, "primn"]), 0)
  # non-prim scaled: factor = (280 - 100) / 200 = 0.9, urc was 200 after clamping
  expect_equal(as.vector(out[, 2025, "urc"]), 180)
  expect_equal(as.vector(out[, 2025, c("secdf", "secdn")]), c(0, 0))
})
