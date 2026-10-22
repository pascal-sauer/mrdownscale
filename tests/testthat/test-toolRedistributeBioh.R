test_that("toolRedistributeBioh redistributes within the region", {
  clusters <- c("regone.1", "regone.2", "regone.3", "regtwo.1")
  bioh <- new.magpie(clusters, years = 2025, names = c("bioh.primf", "bioh.secmf"), fill = 0)
  bioh["regone.1", , ] <- c(60, 120)
  bioh["regone.2", , ] <- c(360, 420)
  bioh["regone.3", , ] <- c(15, 25)
  bioh["regtwo.1", , ] <- c(100, 100)
  calculated <- new.magpie(clusters, years = 2025, names = c("bioh.primf", "bioh.secmf"), fill = 0)
  calculated["regone.2", , ] <- c(200, 400)
  calculated["regone.3", , ] <- c(100, 200)
  calculated["regtwo.1", , ] <- c(100, 100)
  getSets(bioh) <- c("region", "id", "year", "category", "data")
  getSets(calculated) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(
    redistributed <- toolRedistributeBioh(bioh, calculated)
  )
  expect_false(any(grepl("Check failed", conditions$conditions)))
  expect_true(any(grepl("180 kg C yr-1 \\(15%", conditions$conditions)))
  expect_true(any(grepl("from 1 of 4 cluster years \\(25%\\)", conditions$conditions)))

  # the stranded bioh of regone.1 fits into its own region, so it is split
  # proportional to the calculated bioh of its region partners and regtwo stays
  # untouched
  expect_equal(as.vector(redistributed["regone.1", , ]), c(0, 0))
  expect_equal(as.vector(redistributed["regone.2", , ]), c(360 + 40, 420 + 80))
  expect_equal(as.vector(redistributed["regone.3", , ]), c(15 + 20, 25 + 40))
  expect_equal(as.vector(redistributed["regtwo.1", , ]), c(100, 100))
  expect_equal(as.numeric(dimSums(redistributed, c(1, 2, 3))), as.numeric(dimSums(bioh, c(1, 2, 3))))
})

test_that("toolRedistributeBioh falls back to the global pool and drops what cannot move", {
  clusters <- c("regone.1", "regtwo.1")
  bioh <- new.magpie(clusters, years = 2025, names = c("bioh.primf", "bioh.secmf"), fill = 0)
  bioh["regone.1", , ] <- c(60, 120)
  bioh["regtwo.1", , ] <- c(300, 400)
  calculated <- new.magpie(clusters, years = 2025, names = c("bioh.primf", "bioh.secmf"), fill = 0)
  calculated["regtwo.1", , ] <- c(100, 200)
  getSets(bioh) <- c("region", "id", "year", "category", "data")
  getSets(calculated) <- c("region", "id", "year", "category", "data")

  conditions <- captureConditions(
    redistributed <- toolRedistributeBioh(bioh, calculated)
  )
  # regone has no receiver of its own, so the bioh moves over the global
  # pool to regtwo, the only receiver
  expect_false(any(grepl("Check failed", conditions$conditions)))
  expect_true(any(grepl("from 1 of 2 cluster years \\(50%\\)", conditions$conditions)))
  expect_equal(as.vector(redistributed["regone.1", , ]), c(0, 0))
  expect_equal(as.vector(redistributed["regtwo.1", , ]), c(360, 520))

  # if nowhere has receivers, the bioh is dropped and reported as a warning,
  # so the zero bioh on zero harvest area invariant is kept
  conditions <- captureConditions(
    result <- toolRedistributeBioh(bioh, calculated * 0)
  )
  expect_true(any(grepl("dropped stranded bioh without receivers: 880 kg C yr-1 \\(100%", conditions$conditions)))
  expect_true(any(grepl("in 2 of 2 cluster years \\(100%\\)", conditions$conditions)))
  expect_equal(as.numeric(result), rep(0, 4))
})
