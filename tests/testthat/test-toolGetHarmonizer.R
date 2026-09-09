test_that("toolGetHarmonizer returns the requested harmonizer", {
  expect_true(is.function(toolGetHarmonizer("offset")))
  expect_true(is.function(toolGetHarmonizer("fade")))
  expect_true(is.function(toolGetHarmonizer("fadeForest")))
  expect_true(is.function(toolGetHarmonizer("absoluteChanges")))
  expect_error(toolGetHarmonizer("nonexistent"), "harmonizerName")
})
