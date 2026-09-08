test_that("toolGetHarmonizer returns known harmonizers", {
  expect_true(is.function(toolGetHarmonizer("offset")))
  expect_true(is.function(toolGetHarmonizer("fade")))
  expect_true(is.function(toolGetHarmonizer("fadeForest")))
  expect_true(is.function(toolGetHarmonizer("absoluteChanges")))
})

test_that("toolGetHarmonizer errors for an unknown name", {
  expect_error(toolGetHarmonizer("doesNotExist"))
})
