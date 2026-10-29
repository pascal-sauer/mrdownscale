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

test_that("toolParallelHarvestCoverTotals keeps cover targets parallel and allocates onto harvestable land", {
  # regone grows its forest harvest in the input while land stays flat, regtwo
  # is inert: every region's cover target is its target total plus its own
  # input delta, allocated onto the clusters proportional to their harvestable land
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 100, primn = 40, secdf = 50, secdn = 20, pltns = 10),
    regone.2 = list(primf = 40, primn = 20, secdf = 40, secdn = 10, pltns = 10),
    regtwo.1 = list(primf = 30, primn = 10, secdf = 30, secdn = 15, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor, regone.2 = anchor, regtwo.1 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = c(2, 8, 14), pltns = 1, secnf = 1),
                             regone.2 = anchor, regtwo.1 = anchor))
  primHarv <- new.magpie(c("regone.1", "regone.2", "regtwo.1"), years = c(2025, 2030),
                         names = c("primf", "primn"), fill = 0)

  conditions <- captureConditions(
    result <- toolParallelHarvestCoverTotals(xInput, xTarget, primHarv, land, 2020)
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))

  # forest target 6 + (12 - 6) = 12 in 2025 for regone and 3 for regtwo,
  # other land stays at its anchors of 2 and 1, the global totals balanced
  expect_equal(as.vector(result$coverTarget$forest["regone", 2025, ]), 12)
  expect_equal(as.vector(result$coverTarget$other["regtwo", 2030, ]), 1)
  expect_equal(as.vector(result$coverGlobalBalanced[, 2025, ]), 15 + 3)
  expect_identical(result$covers,
                   list(forest = list(prim = "primf", secd = c("secdf", "pltns")),
                        other = list(prim = "primn", secd = "secdn")))

  # the allocation follows the harvestable land caps (60 : 50 within regone),
  # the cluster shares are their item caps of the cover caps
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, "secdf"]), 12 * 50 / 110, tolerance = 10^-8)
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, "pltns"]), 12 * 10 / 110, tolerance = 10^-8)
  expect_equal(as.vector(result$secdHarv["regone.2", 2025, "secdf"]), 12 * 40 / 110, tolerance = 10^-8)
  expect_equal(as.vector(result$secdHarv["regtwo.1", 2025, "secdn"]), 1)
  expect_equal(max(result$excessSecdHarvest), 0)
})

test_that("toolParallelHarvestCoverTotals reports the excess that not even the other cover can host", {
  # regone grows its forest harvest to 23 while its forest land caps at 12 and
  # the other land cover takes just 1 of the 11 excess: what remains unhostable
  # is dropped at the caps and returned as an excess of 10
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
  primHarv <- new.magpie(c("regone.1", "regone.2"), years = c(2025, 2030),
                         names = c("primf", "primn"), fill = 0)

  conditions <- captureConditions(
    result <- toolParallelHarvestCoverTotals(xInput, xTarget, primHarv, land, 2020)
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))

  expect_equal(as.numeric(dimSums(result$excessSecdHarvest, 1)[, 2025, ]), 10)
  expect_equal(as.numeric(dimSums(result$excessSecdHarvest, 1)[, 2030, ]), 10)

  # both clusters end exactly on their caps, the hostable other land target
  # fills its clusters to their caps too
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(result$secdHarv["regone.2", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, "secdn"]), 1)
  expect_equal(as.vector(result$coverTarget$forest["regone", 2025, ]), 23)
  expect_equal(as.vector(result$coverGlobalBalanced[, 2025, ]), 23 + 1)
})
