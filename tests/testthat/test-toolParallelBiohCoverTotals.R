test_that("toolParallelBiohCoverTotals moves stranded bioh over the global pool, failing the regional check", {
  # regone.1 targets 2500 + (2100 - 3500) = 1100 forest bioh in 2025/2030 but
  # its recalculated bioh is zero, so it hosts nothing and its target
  # moves over the global pool to regtwo.1 with its recalculated 1000/1000
  cells <- c("regone.1", "regtwo.1")
  items <- c("bioh.primf", "bioh.secmf")

  xInput <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput["regone.1", 2020, ] <- c(2000, 1500)
  xInput["regone.1", c(2025, 2030), "bioh.primf"] <- 500
  xInput["regone.1", c(2025, 2030), "bioh.secmf"] <- 1600
  xInput["regtwo.1", , ] <- 1000
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget["regone.1", , "bioh.primf"] <- 1000
  xTarget["regone.1", , "bioh.secmf"] <- 1500
  xTarget["regtwo.1", , ] <- 1000
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  biohCalculated <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 1000)
  biohCalculated["regone.1", c(2025, 2030), ] <- 0

  conditions <- captureConditions(
    result <- toolParallelBiohCoverTotals(xInput, xTarget, biohCalculated, 2020)
  )
  expect_true(any(grepl("Check failed.*Regional forest bioh totals are parallel", conditions$conditions)))
  expect_false(any(grepl("Check failed.*Global bioh totals are parallel", conditions$conditions)))

  # regtwo.1 hosts its own 2000 plus the pooled 1100 of regone, split by its
  # recalculated item weights 1000/1000
  expect_equal(as.vector(result["regone.1", c(2025, 2030), ]), rep(0, 4))
  expect_equal(as.vector(result["regtwo.1", c(2025, 2030), ]), rep((2000 + 1100) / 2, 4))

  # the target years up to the harmonization year are returned unchanged
  expect_equal(as.vector(result[, 2010, ]), as.vector(xTarget[, 2010, ]))
  expect_equal(as.vector(result[, 2020, ]), as.vector(xTarget[, 2020, ]))
})

test_that("toolParallelBiohCoverTotals allocates the global cover target first when a regional one is negative", {
  # regone's other land target 100 + (1 - 300) = -199 is floored to zero
  # while regtwo grows to 400: scaling the positive target onto the raw
  # global one 201 keeps the global cover total parallel, the forest cover with
  # positive targets everywhere is untouched by the balancing
  cells <- c("regone.1", "regtwo.1")
  items <- c("bioh.primf", "bioh.secnf")

  xInput <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput[, , "bioh.primf"] <- 1000
  xInput["regone.1", , "bioh.secnf"] <- c(300, 1, 1)
  xInput["regtwo.1", , "bioh.secnf"] <- c(100, 400, 600)
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , "bioh.primf"] <- 1000
  xTarget["regone.1", , "bioh.secnf"] <- 100
  xTarget["regtwo.1", , "bioh.secnf"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  biohCalculated <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 1000)

  conditions <- captureConditions(
    result <- toolParallelBiohCoverTotals(xInput, xTarget, biohCalculated, 2020)
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # regone hosts nothing of the cover whose regional target turned
  # negative, regtwo carries the whole raw global target 201 and 401
  expect_equal(as.vector(result["regone.1", c(2025, 2030), "bioh.secnf"]), c(0, 0))
  expect_equal(as.vector(result["regtwo.1", 2025, "bioh.secnf"]), 201, tolerance = 10^-9)
  expect_equal(as.vector(result["regtwo.1", 2030, "bioh.secnf"]), 401, tolerance = 10^-9)

  # the forest cover stays exactly parallel
  expect_equal(as.vector(result["regone.1", 2025, "bioh.primf"]), 1000, tolerance = 10^-9)
  expect_equal(as.vector(result["regtwo.1", 2030, "bioh.primf"]), 1000, tolerance = 10^-9)
})

test_that("toolParallelBiohCoverTotals charges a negative global cover target to the other covers", {
  # both regions' other land targets are negative, so the global other
  # land target -298 turns negative and is floored to zero, charged to the
  # forest cover whose regions are uniformly scaled by 1702 / 2000
  cells <- c("regone.1", "regtwo.1")
  items <- c("bioh.primf", "bioh.secnf")

  xInput <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 0)
  xInput[, , "bioh.primf"] <- 1000
  xInput["regone.1", , "bioh.secnf"] <- c(300, 1, 1)
  xInput["regtwo.1", , "bioh.secnf"] <- c(200, 1, 1)
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  xTarget <- new.magpie(cells, years = c(2010, 2020), names = items, fill = 0)
  xTarget[, , "bioh.primf"] <- 1000
  xTarget["regone.1", , "bioh.secnf"] <- 100
  xTarget["regtwo.1", , "bioh.secnf"] <- 100
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  biohCalculated <- new.magpie(cells, years = c(2020, 2025, 2030), names = items, fill = 1000)

  conditions <- captureConditions(
    result <- toolParallelBiohCoverTotals(xInput, xTarget, biohCalculated, 2020)
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))

  expect_equal(as.vector(result[, c(2025, 2030), "bioh.secnf"]), rep(0, 4))
  # the forest cover absorbs the charged loss: 1000 * 1702 / 2000 per region
  expect_equal(as.vector(result["regone.1", 2025, "bioh.primf"]), 1000 * (1702 / 2000), tolerance = 10^-6)
  expect_equal(as.vector(result["regtwo.1", 2030, "bioh.primf"]), 1000 * (1702 / 2000), tolerance = 10^-6)
  # the global total over all covers stays parallel to the input
  expect_equal(as.numeric(dimSums(result[, 2025, items], c(1, 2, 3))), 1702, tolerance = 10^-6)
  expect_equal(as.numeric(dimSums(result[, 2030, items], c(1, 2, 3))), 1702, tolerance = 10^-6)
})
