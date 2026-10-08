test_that("toolHandleNegatives rescales cells to their target area", {
  x <- new.magpie(1:2, names = c("a", "b", "c"), fill = c(0.5, -0.2, 0.7, 1, 1, 1))
  expect_silent(res <- toolHandleNegatives(x))
  expect_equal(dimSums(res, 3), dimSums(x, 3))
  expect_true(all(res >= 0))
  expect_equal(as.numeric(toolHandleNegatives(x[, , "c"])), as.numeric(x[, , "c"]))
  # cells whose area was fully cancelled by negatives stay zero
  x2 <- new.magpie(1:2, names = c("a", "b"), fill = c(0, 0.5, 0, -0.5))
  expect_equal(as.numeric(toolHandleNegatives(x2)), rep(0, 4))
  # every year keeps its own total area
  x3 <- new.magpie(1, years = c(1995, 2000), names = c("a", "b"), fill = c(0.6, 1.2, 0.4, -0.4))
  expect_equal(as.numeric(dimSums(toolHandleNegatives(x3), 3)), c(1, 0.8))
  # an explicit target area is recycled to all years
  target <- new.magpie(1, names = "target", fill = 0.9)
  expect_equal(as.numeric(dimSums(toolHandleNegatives(x3, targetArea = target), 3)), rep(0.9, 2))
  # a negative target area is clipped to zero, zeroing the affected cells
  x4 <- new.magpie(1, names = c("a", "b"), fill = c(0.5, -1))
  expect_equal(as.numeric(toolHandleNegatives(x4, targetArea = pmax(dimSums(x4, 3), 0))), rep(0, 2))
  # only years with a negative target area are zeroed
  x5 <- new.magpie(1, years = c(1995, 2000), names = c("a", "b"), fill = c(0.5, 0.5, 0, -1))
  expect_equal(as.numeric(toolHandleNegatives(x5, targetArea = pmax(dimSums(x5, 3), 0))), c(0.5, 0, 0, 0))
})

test_that("toolHandleNegatives scales lastScaled categories only as a last resort", {
  x <- new.magpie(1, names = c("a", "b", "c"), fill = c(1, 1, 0.5))
  target <- new.magpie(1, names = "target", fill = 0.75)
  # the excess of 1.75 is fully absorbed by a and b, c stays untouched
  expect_equal(as.numeric(toolHandleNegatives(x, targetArea = target, lastScaled = "c")),
               c(0.125, 0.125, 0.5))
  # a and b are scaled down to zero first, c only absorbs the remaining 0.25
  expect_equal(as.numeric(toolHandleNegatives(x, targetArea = new.magpie(1, names = "target", fill = 0.25),
                                              lastScaled = "c")),
               c(0, 0, 0.25))
  # without lastScaled the old proportional behavior is unchanged
  expect_equal(as.numeric(toolHandleNegatives(x, targetArea = target)), c(0.3, 0.3, 0.15))
})

test_that("toolHandleNegatives errors on invalid input", {
  x <- new.magpie(1, names = c("a", "b"), fill = c(0.5, -1))
  expect_error(toolHandleNegatives(x), "all(targetArea >= 0) is not TRUE", fixed = TRUE)
  expect_error(toolHandleNegatives(x, targetArea = new.magpie(1, names = "target", fill = -5)),
               "all(targetArea >= 0) is not TRUE", fixed = TRUE)
  # negative target areas are clipped by the caller, not handled here
  expect_error(toolHandleNegatives(x, allowNegativeTarget = TRUE), "unused argument")
  # scaling up is not allowed
  x2 <- new.magpie(1, names = c("a", "b"), fill = c(1, 1))
  expect_error(toolHandleNegatives(x2, targetArea = new.magpie(1, names = "target", fill = 3)),
               "all(fact <= 1 + tolerance) is not TRUE", fixed = TRUE)
  # missing values are rejected, also in later years
  xna <- new.magpie(1, names = c("a", "b"), fill = c(0.5, 1))
  xna[1, 1, 1] <- NA
  expect_error(toolHandleNegatives(xna), "all(targetArea >= 0) is not TRUE", fixed = TRUE)
  xna2 <- new.magpie(1, years = c(1995, 2000), names = c("a", "b"), fill = c(0.6, 1, 0.4, 1))
  xna2[1, 2, 1] <- NA
  expect_error(toolHandleNegatives(xna2), "all(targetArea >= 0) is not TRUE", fixed = TRUE)
})
