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

  # the allocation follows the harvestable land caps (60 : 50 within regone),
  # the cluster shares are their item caps of the cover caps, the forest
  # target being 6 + (12 - 6) = 12 in 2025 for regone and the other land
  # staying at its anchor of 2
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
  # the self checks report the parallelism broken by the caps on all margins
  expect_true(any(grepl("Check failed.*Regional forest wood harvest area is parallel", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Regional other wood harvest area is parallel", conditions$conditions)))
  expect_true(any(grepl("Check failed.*Global wood harvest area is parallel", conditions$conditions)))

  expect_equal(as.numeric(dimSums(result$excessSecdHarvest, 1)[, 2025, ]), 10)
  expect_equal(as.numeric(dimSums(result$excessSecdHarvest, 1)[, 2030, ]), 10)

  # both clusters end exactly on their caps, the hostable other land target
  # fills its clusters to their caps too
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(result$secdHarv["regone.2", 2025, c("secdf", "pltns")]), c(5, 1))
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, "secdn"]), 1)
})

test_that("toolParallelHarvestCoverTotals places a stranded target into another region's room", {
  # regtwo targets 3 forest harvest but has no forest land at all, while
  # regone's own target 4 leaves 6 room on its cluster: the whole stranded 3
  # fits into that room, so nothing is lost and no excess is reported
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 0, secdf = 10, secdn = 0, pltns = 0),
    regtwo.1 = list(primf = 20, primn = 0, secdf = 0, secdn = 0, pltns = 0)
  ))

  xTarget <- makeHarvest(c(2010, 2020), list(
    regone.1 = list(primf = 0, primn = 0, secyf = 4, secmf = 0, pltns = 0, secnf = 0),
    regtwo.1 = list(primf = 0, primn = 0, secyf = 3, secmf = 0, pltns = 0, secnf = 0)
  ))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = 4, pltns = 0, secnf = 0),
                             regtwo.1 = list(secyf = 3, pltns = 0, secnf = 0)))
  primHarv <- new.magpie(c("regone.1", "regtwo.1"), years = c(2025, 2030),
                         names = c("primf", "primn"), fill = 0)

  conditions <- captureConditions(
    result <- toolParallelHarvestCoverTotals(xInput, xTarget, primHarv, land, 2020)
  )
  expect_true(any(grepl("Check failed.*harvestable land left to hold it", conditions$conditions)))

  # all 7 of the global forest target land on regone.1 within its 10 cap
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, "secdf"]), 7, tolerance = 10^-8)
  expect_equal(as.vector(result$secdHarv["regtwo.1", 2025, "secdf"]), 0)
  expect_equal(max(result$excessSecdHarvest), 0)
})

test_that("toolParallelHarvestCoverTotals reports the stranded target that no room can host", {
  # both regions target forest harvest but there is no forest land at all: the
  # whole stranded forest target finds no room anywhere and must be reported
  # as excess spread over the clusters instead of vanishing silently
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 20, primn = 0, secdf = 0, secdn = 0, pltns = 0),
    regtwo.1 = list(primf = 20, primn = 0, secdf = 0, secdn = 0, pltns = 0)
  ))

  xTarget <- makeHarvest(c(2010, 2020), list(
    regone.1 = list(primf = 0, primn = 0, secyf = 4, secmf = 0, pltns = 0, secnf = 0),
    regtwo.1 = list(primf = 0, primn = 0, secyf = 3, secmf = 0, pltns = 0, secnf = 0)
  ))
  xInput <- makeHarvest(c(2020, 2025, 2030),
                        list(regone.1 = list(secyf = 4, pltns = 0, secnf = 0),
                             regtwo.1 = list(secyf = 3, pltns = 0, secnf = 0)))
  primHarv <- new.magpie(c("regone.1", "regtwo.1"), years = c(2025, 2030),
                         names = c("primf", "primn"), fill = 0)

  conditions <- captureConditions(
    result <- toolParallelHarvestCoverTotals(xInput, xTarget, primHarv, land, 2020)
  )
  expect_true(any(grepl("Check failed.*harvestable land left to hold it", conditions$conditions)))

  # nothing can be allocated, the full target of 7 is reported as excess
  expect_equal(as.vector(result$secdHarv["regone.1", 2025, "secdf"]), 0)
  expect_equal(as.numeric(dimSums(result$excessSecdHarvest, 1)[, 2025, ]), 7)
  expect_equal(as.numeric(dimSums(result$excessSecdHarvest, 1)[, 2030, ]), 7)
})

test_that("toolParallelHarvestCoverTotals checks its own allocation for parallelism", {
  # the land-derived primary forest harvest of 10 alone outweighs the forest
  # cover target of 3: the tool itself must report the broken regional and
  # global parallelism of its allocation
  years <- c(2010, 2020, 2025, 2030)
  land <- makeLand(years, list(
    regone.1 = list(primf = 100, primn = 40, secdf = 50, secdn = 20, pltns = 10)
  ))

  anchor <- list(primf = 0, primn = 0, secyf = 2, secmf = 0, pltns = 1, secnf = 1)
  xTarget <- makeHarvest(c(2010, 2020), list(regone.1 = anchor))
  xInput <- makeHarvest(c(2020, 2025, 2030), list(regone.1 = anchor))
  primHarv <- new.magpie("regone.1", years = c(2025, 2030), names = c("primf", "primn"), fill = 0)
  primHarv[, , "primf"] <- 10

  conditions <- captureConditions(
    result <- toolParallelHarvestCoverTotals(xInput, xTarget, primHarv, land, 2020)
  )
  expect_true(any(grepl("Check failed.*Regional forest wood harvest area is parallel to input",
                        conditions$conditions)))
  expect_true(any(grepl("Check failed.*Global wood harvest area is parallel to input",
                        conditions$conditions)))
  expect_identical(names(result), c("secdHarv", "excessSecdHarvest"))
})
