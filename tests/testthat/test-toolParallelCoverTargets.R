makeCover <- function(regionValues, years = c(2020, 2025, 2030)) {
  x <- new.magpie(names(regionValues), years = years, names = "cover", fill = 0)
  for (region in names(regionValues)) {
    x[region, , ] <- regionValues[[region]]
  }
  return(x)
}

test_that("toolParallelCoverTargets keeps positive cover targets untouched and sums them globally", {
  # every regional target is anchor plus input delta and stays positive, so
  # neither the flooring nor the balancing changes it
  inputCover <- list(forest = makeCover(list(regone = c(100, 120, 140), regtwo = c(50, 50, 50))),
                     other = makeCover(list(regone = c(10, 12, 15), regtwo = c(20, 20, 20))))
  targetCover <- list(forest = makeCover(list(regone = 100, regtwo = 50)),
                      other = makeCover(list(regone = 10, regtwo = 20)))

  result <- toolParallelCoverTargets(inputCover, targetCover, 2020)

  expect_equal(as.vector(result$coverTargets$forest["regone", c(2025, 2030), ]), c(120, 140))
  expect_equal(as.vector(result$coverTargets$forest["regtwo", c(2025, 2030), ]), c(50, 50))
  expect_equal(as.vector(result$coverTargets$other["regone", c(2025, 2030), ]), c(12, 15))
  expect_equal(as.vector(result$coverTargets$other["regtwo", c(2025, 2030), ]), c(20, 20))
  expect_equal(as.numeric(result$coverGlobalBalanced[, 2025, ]), 202)
  expect_equal(as.numeric(result$coverGlobalBalanced[, 2030, ]), 225)
})

test_that("toolParallelCoverTargets floors negative regional targets and keeps the global cover total parallel", {
  # regone's forest target 100 + (1 - 300) = -199 is floored to zero while
  # regtwo's 400 is scaled onto the raw global target 201, keeping the global
  # forest cover total input parallel
  inputCover <- list(forest = makeCover(list(regone = c(300, 1), regtwo = c(100, 400)),
                                        years = c(2020, 2025)),
                     other = makeCover(list(regone = c(10, 10), regtwo = c(20, 20)),
                                       years = c(2020, 2025)))
  targetCover <- list(forest = makeCover(list(regone = 100, regtwo = 100), years = c(2020, 2025)),
                      other = makeCover(list(regone = 10, regtwo = 20), years = c(2020, 2025)))

  result <- toolParallelCoverTargets(inputCover, targetCover, 2020)

  expect_equal(as.vector(result$coverTargets$forest["regone", 2025, ]), 0)
  expect_equal(as.vector(result$coverTargets$forest["regtwo", 2025, ]), 201, tolerance = 10^-9)
  expect_equal(as.numeric(dimSums(result$coverTargets$forest, 1)[, 2025, ]), 201, tolerance = 10^-9)
  # the other cover stays untouched, the sum over covers stays input parallel
  expect_equal(as.vector(result$coverTargets$other["regone", 2025, ]), 10)
  expect_equal(as.numeric(dimSums(result$coverTargets$forest + result$coverTargets$other,
                                  1)[, 2025, ]), 231, tolerance = 10^-9)
  expect_equal(as.numeric(result$coverGlobalBalanced[, 2025, ]), 231, tolerance = 10^-9)
})

test_that("toolParallelCoverTargets floors a negative global cover target at zero and charges the other covers", {
  # both regions' other land targets are negative (-199 and 100 + (1 - 200)
  # = -99), so the global other target -298 is floored to zero and its loss
  # charged to the forest cover whose regions are uniformly scaled by
  # 1702 / 2000, keeping the global total over all covers input parallel
  inputCover <- list(forest = makeCover(list(regone = c(1000, 1000), regtwo = c(1000, 1000)),
                                        years = c(2020, 2025)),
                     other = makeCover(list(regone = c(300, 1), regtwo = c(200, 1)),
                                       years = c(2020, 2025)))
  targetCover <- list(forest = makeCover(list(regone = 1000, regtwo = 1000), years = c(2020, 2025)),
                      other = makeCover(list(regone = 100, regtwo = 100), years = c(2020, 2025)))

  result <- toolParallelCoverTargets(inputCover, targetCover, 2020)

  expect_equal(as.vector(result$coverTargets$other["regone", 2025, ]), 0)
  expect_equal(as.vector(result$coverTargets$other["regtwo", 2025, ]), 0)
  expect_equal(as.vector(result$coverTargets$forest["regone", 2025, ]), 1000 * (1702 / 2000),
               tolerance = 10^-9)
  expect_equal(as.vector(result$coverTargets$forest["regtwo", 2025, ]), 1000 * (1702 / 2000),
               tolerance = 10^-9)
  expect_equal(as.numeric(dimSums(result$coverTargets$forest + result$coverTargets$other,
                                  1)[, 2025, ]), 1702, tolerance = 10^-9)
  expect_equal(as.numeric(result$coverGlobalBalanced[, 2025, ]), 1702, tolerance = 10^-9)
})

test_that("toolParallelCoverTargets handles all-zero targets without producing NaN", {
  # the 0 / 0 scale and balance factor are replaced, no NaN leaks out
  inputCover <- list(forest = makeCover(list(regone = c(0, 0), regtwo = c(0, 0)),
                                        years = c(2020, 2025)),
                     other = makeCover(list(regone = c(0, 0), regtwo = c(0, 0)),
                                       years = c(2020, 2025)))
  targetCover <- list(forest = makeCover(list(regone = 0, regtwo = 0), years = c(2020, 2025)),
                      other = makeCover(list(regone = 0, regtwo = 0), years = c(2020, 2025)))

  result <- toolParallelCoverTargets(inputCover, targetCover, 2020)

  expect_equal(as.vector(result$coverTargets$forest), rep(0, 2))
  expect_equal(as.vector(result$coverTargets$other), rep(0, 2))
  expect_equal(as.numeric(result$coverGlobalBalanced), 0)
})
