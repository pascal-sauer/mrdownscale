# helper for toolReportAreaDeviation, which uses messages
captureConditions <- function(expr) {
  conditions <- character()
  value <- withCallingHandlers(
    expr,
    message = function(m) {
      conditions <<- unique(c(conditions, conditionMessage(m)))
      invokeRestart("muffleMessage")
    },
    warning = function(w) {
      conditions <<- unique(c(conditions, conditionMessage(w)))
      invokeRestart("muffleWarning")
    }
  )
  return(list(value = value, conditions = conditions))
}

test_that("toolReportAreaDeviation reports positive and negative deviations separately", {
  xOut <- new.magpie(c("reg.a", "reg.b"), years = c(2020, 2025, 2030),
                     names = c("c3ann_rainfed", "urban"), fill = 0)
  xOut["reg.a", 2020, "c3ann_rainfed"] <- 10
  xOut["reg.b", 2020, "c3ann_rainfed"] <- 10
  xOut["reg.a", 2025, "c3ann_rainfed"] <- 12
  xOut["reg.b", 2025, "c3ann_rainfed"] <- 7
  xOut["reg.a", 2030, "c3ann_rainfed"] <- 11
  xOut["reg.b", 2030, "c3ann_rainfed"] <- 9

  xRaw <- new.magpie(c("reg.a", "reg.b"), years = c(2025, 2030),
                     names = c("c3ann_rainfed", "urban"), fill = 0)
  xRaw["reg.a", , "c3ann_rainfed"] <- c(14, 13)
  xRaw["reg.b", , "c3ann_rainfed"] <- c(5, 8)

  run <- captureConditions(toolReportAreaDeviation(
    xRaw, xOut,
    groups = list("forest and other land" = character(0), cropland = "c3ann_rainfed",
                  "pasture and rangeland" = character(0))
  ))

  cropGroup <- gsub("[[:space:]]", "", grep("group \"cropland\"", run$conditions, value = TRUE))
  expect_length(cropGroup, 1)
  expect_length(grep("forest", run$conditions), 0)
  # deviations of -2, -2, +2, +1 Mha accumulate to -4/+3 Mha over the two
  # timesteps, means -2/+1.5 Mha, 39 Mha of output area give -10.3/+7.69 %,
  # net change -1 Mha is -0.5 Mha per timestep and -2.56 %
  expect_match(cropGroup, "c3ann_rainfed,-2\\.56,7\\.69,-10\\.3,-0\\.5,1\\.5,-2")
  expect_match(cropGroup, "total,-2\\.56,7\\.69,-10\\.3,-0\\.5,1\\.5,-2")
})

test_that("toolReportAreaDeviation sorts by gross deviation and sums group totals", {
  xOut <- new.magpie("reg.a", years = c(2020, 2025, 2030),
                     names = c("c3ann_rainfed", "c4ann_rainfed", "urban"), fill = 0)
  xOut["reg.a", , "c3ann_rainfed"] <- c(10, 10, 10)
  xOut["reg.a", , "c4ann_rainfed"] <- c(20, 17, 17)

  xRaw <- new.magpie("reg.a", years = c(2025, 2030),
                     names = c("c3ann_rainfed", "c4ann_rainfed", "urban"), fill = 0)
  # c4ann_rainfed deviates by +7 Mha gross, more than the +2 of c3ann_rainfed,
  # so it must be listed first although it comes second in the group
  xRaw["reg.a", , "c3ann_rainfed"] <- c(11, 9)
  xRaw["reg.a", , "c4ann_rainfed"] <- c(21, 20)

  run <- captureConditions(toolReportAreaDeviation(
    xRaw, xOut,
    groups = list("forest and other land" = character(0),
                  cropland = c("c3ann_rainfed", "c4ann_rainfed"),
                  "pasture and rangeland" = character(0))
  ))

  cropGroup <- gsub("[[:space:]]", "", grep("group \"cropland\"", run$conditions, value = TRUE))
  expect_length(cropGroup, 1)
  # c3ann_rainfed: net 0, +1/-1 over output area 20; c4ann_rainfed: net -7,
  # 0/-7 over output area 34; totals: net -7, +1/-8 over output area 54
  expect_match(cropGroup, "c4ann_rainfed,-20\\.6,0,-20\\.6,-3\\.5,0,-3\\.5")
  expect_match(cropGroup, "c3ann_rainfed,0,5,-5,0,0\\.5,-0\\.5")
  expect_match(cropGroup, "total,-13,1\\.85,-14\\.8,-3\\.5,0\\.5,-4")
  expect_true(gregexpr("c4ann_rainfed", cropGroup)[[1]][1] <
                gregexpr("c3ann_rainfed", cropGroup)[[1]][1])
})

test_that("toolReportAreaDeviation rejects item order mismatch", {
  xOut <- new.magpie("reg.a", years = c(2020, 2025),
                     names = c("c3ann_rainfed", "urban"), fill = 1)
  # magclass arithmetic aligns silently and keeps the item order of the first
  # operand, so without an assertion the values would be relabeled wrongly
  xRaw <- new.magpie("reg.a", years = 2025, names = c("urban", "c3ann_rainfed"), fill = 1)
  expect_error(toolReportAreaDeviation(xRaw, xOut, groups = list(cropland = "c3ann_rainfed")),
               "identical\\(getItems")
})

test_that("toolReportAreaDeviation handles zero output area", {
  xOut <- new.magpie("reg.a", years = c(2020, 2025, 2030),
                     names = c("c4ann_rainfed", "urban"), fill = 0)
  xRaw <- new.magpie("reg.a", years = c(2025, 2030),
                     names = c("c4ann_rainfed", "urban"), fill = 0)
  # the output has no area of this variable, so only the net negative
  # deviation is infinite and the absent positive one is n/a, not NaN
  xRaw["reg.a", , "c4ann_rainfed"] <- c(1, 0)

  run <- captureConditions(toolReportAreaDeviation(
    xRaw, xOut,
    groups = list("forest and other land" = character(0), cropland = "c4ann_rainfed",
                  "pasture and rangeland" = character(0))
  ))

  cropGroup <- gsub("[[:space:]]", "", grep("group \"cropland\"", run$conditions, value = TRUE))
  expect_length(cropGroup, 1)
  expect_match(cropGroup, "c4ann_rainfed,-Inf,n/a,-Inf,-0\\.5,0,-0\\.5")
})
