test_that("toolProjectRegionalDelta projects regional deltas proportional to the share weight", {
  clusters <- c("regone.1", "regone.2", "regtwo.1", "regthree.1")

  xTarget <- new.magpie(clusters, years = c(2010, 2020), names = "amount.c3ann", fill = 0)
  xTarget["regone.1", , ] <- 2
  xTarget["regone.2", , ] <- 2
  xTarget["regtwo.1", , ] <- 0.5
  xTarget["regthree.1", , ] <- 3
  getSets(xTarget) <- c("region", "id", "year", "category", "data")

  xInput <- new.magpie(clusters, years = c(2020, 2025), names = "amount.c3ann", fill = 0)
  xInput["regone.1", , ] <- c(10, 14)
  xInput["regone.2", , ] <- c(1, 1)
  xInput["regtwo.1", , ] <- c(1, 0)
  xInput["regthree.1", , ] <- c(1, 5)
  getSets(xInput) <- c("region", "id", "year", "category", "data")

  shareWeight <- new.magpie(clusters, years = 2025, names = "c3ann", fill = 0)
  shareWeight["regone.1", , ] <- 0.003
  shareWeight["regone.2", , ] <- 20
  shareWeight["regtwo.1", , ] <- 100

  projected <- toolProjectRegionalDelta(xTarget, xInput, shareWeight, 2020)

  expect_equal(getYears(projected, as.integer = TRUE), 2025)
  expect_true(setequal(getItems(projected, 3), "amount.c3ann"))

  # the +4 regional delta of regone is distributed onto its clusters proportional
  # to their share weights, on top of their anchors in the harmonization year
  expect_equal(as.vector(projected["regone.1", , ]), 2 + 4 * (0.003 / 20.003), tolerance = 10^-9)
  expect_equal(as.vector(projected["regone.2", , ]), 2 + 4 * (20 / 20.003), tolerance = 10^-9)
  # projections below zero are clamped
  expect_equal(as.vector(projected["regtwo.1", , ]), 0)
  # a region without any share weight keeps its anchor despite its delta
  expect_equal(as.vector(projected["regthree.1", , ]), 3)

  # the regional total of anchor plus delta is preserved where nothing was clamped
  expect_equal(as.numeric(dimSums(projected["regone", , ], c(1, 3))), 4 + 4)
})
