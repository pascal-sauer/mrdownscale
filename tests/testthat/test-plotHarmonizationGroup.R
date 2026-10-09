test_that("plotHarmonizationGroup writes a png and returns its path", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(Year = rep(2019:2021, 4),
                  Item = rep(c("a", "b"), each = 3, times = 2),
                  Value = seq_len(12),
                  Method = rep(c("absoluteChanges", "inputUnharmonized"), each = 6))
  dir <- withr::local_tempdir()
  file <- mrdownscale:::plotHarmonizationGroup("land", "test", d, ylab = "x",
                                               harmonizationPeriod = c(2020, 2050), plotdir = dir)
  expect_identical(file, file.path(dir, "plot-harmonized-land-test.png"))
  expect_true(file.exists(file))
})

test_that("plotHarmonizationGroup works with a fade method (second vline)", {
  skip_if_not_installed("ggplot2")
  d <- data.frame(Year = rep(2019:2021, 4),
                  Item = rep(c("a", "b"), each = 3, times = 2),
                  Value = seq_len(12),
                  Method = rep(c("fadeForest", "fade"), each = 6))
  dir <- withr::local_tempdir()
  file <- mrdownscale:::plotHarmonizationGroup("nonland", "test", d, ylab = "x",
                                               harmonizationPeriod = c(2020, 2050), plotdir = dir)
  expect_identical(file, file.path(dir, "plot-harmonized-nonland-test.png"))
  expect_true(file.exists(file))
})
