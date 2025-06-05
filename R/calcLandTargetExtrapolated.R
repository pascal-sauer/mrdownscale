#' calcLandTargetExtrapolated
#'
#' Aggregated low resolution target data is extrapolated to the given years
#' using toolExtrapolate and normalized afterwards, so that the total sum over
#' all land types is unchanged.
#' To account for the relationship between wood
#' harvest area and primary land (which is no longer primary once it has been
#' harvested) wood harvest area is calculated here even though it is otherwise
#' considered a nonland variable. The share of woody land that was
#' harvested in the historical period is calculated and then multiplied
#' by the maximum possible harvest in the extrapolation period. Primary land
#' is then converted to secondary land so the total reduction equals the area
#' that was harvested.
#'
#' @param input character, name of the input data set
#' @param target character, name of the target data set
#' @param harmonizationPeriod Two integer values, will extrapolate to all years
#' present in input data between harmonization start and end year
#' @return extrapolated land target data, if calcOutput is called with
#' supplementary = TRUE and target is luh2mod wood harvest area is also returned
#' @author Pascal Sauer
calcLandTargetExtrapolated <- function(input, target, harmonizationPeriod) {
  xInput <- calcOutput("LandInputRecategorized", input = input, target = target, aggregate = FALSE)
  inputYears <- getYears(xInput, as.integer = TRUE)
  transitionYears <- inputYears[inputYears > harmonizationPeriod[1] & inputYears < harmonizationPeriod[2]]

  xTarget <- calcOutput("LandTargetLowRes", input = input, target = target, aggregate = FALSE)
  histYears <- getYears(xTarget, as.integer = TRUE)

  # ---------- extrapolate -------------
  exTarget <- toolExtrapolate(xTarget, transitionYears)
  toolWarnIfExpansion(exTarget, c("primf", "primn"))
  exTarget[exTarget < 0] <- 0

  # ---------- normalize -------------
  # normalize exTarget so that its total sum over all layers agrees for all time steps
  # with the sum over all layers in target in the harmonization year (e.g. makes sure
  # that the land changes in a land data set do not alter the total sum of land.)
  targetArea <- dimSums(setYears(xTarget[, max(histYears), ], NULL), dim = 3)
  exTarget <- exTarget * targetArea / dimSums(exTarget, dim = 3)
  exTarget[is.na(exTarget)] <- 0
  out <- mbind(xTarget, exTarget)

  # normalization can introduce prim expansion
  out <- toolReplaceExpansion(out, "primf", "secdf", noteThreshold = 100, warnThreshold = 100)
  out <- toolReplaceExpansion(out, "primn", "secdn", noteThreshold = 100, warnThreshold = 100)

  harvest <- NULL
  if (target %in% c("luh2mod", "luh3")) {
    # ------- calculate wood harvest shares -------
    harvestHist <- calcOutput("NonlandTargetLowRes", input = input, target = target, aggregate = FALSE)
    harvestHist <- harvestHist[, , endsWith(getItems(harvestHist, 3), "wood_harvest_area")]
    maxHarvestHist <- toolMaxHarvestPerYear(xTarget)
    stopifnot(setequal(getItems(harvestHist, 1), getItems(maxHarvestHist, 1)),
              setequal(getItems(harvestHist[, -1, ], 2), getItems(maxHarvestHist, 2)),
              setequal(getItems(harvestHist, 3), getItems(maxHarvestHist, 3)))

    # calculate share: wood harvest area / max possible harvest
    harvestShare <- dimSums(harvestHist[, -1, ], 2) / dimSums(maxHarvestHist, 2)
    harvestShare[is.na(harvestShare)] <- 0
    harvestShare[harvestShare > 1] <- 1

    secymf <- paste0(c("secyf", "secmf"), "_wood_harvest_area")
    normalization <- dimSums(harvestShare[, , secymf], 3) + 10^-10 # + 10^-10 to ensure sum <= 1
    normalization[normalization < 1] <- 1
    harvestShare[, , secymf] <- harvestShare[, , secymf] / normalization

    stopifnot(0 <= harvestShare, harvestShare <= 1,
              dimSums(harvestShare[, , secymf], 3) <= 1)

    # calculate wood harvest area, reduce primf and primn so they are consistent with harvest
    harvest <- add_columns(harvestHist, paste0("y", transitionYears), dim = 2)

    timestepLength <- unique(diff(transitionYears))
    stopifnot(identical(getYears(harvest), getYears(out)),
              length(timestepLength) == 1, timestepLength > 0)
    for (i in match(transitionYears, getYears(out, as.integer = TRUE))) {
      primfn <- c("primf", "primn")
      secdfn <- c("secdf", "secdn")

      # in toolMaxHarvestPerYear out[, i, ] is only used to determine timestepLength and then thrown away
      harvest[, i, ] <- harvestShare * toolMaxHarvestPerYear(out[, c(i - 1, i), ])

      harvestAgg <- toolAggregateWoodHarvest(harvest[, i, ])
      maxPossiblePrim <- out[, i - 1, primfn] - timestepLength * harvestAgg[, , primfn]

      # more prim than in previous timestep is not possible, convert to secd
      toSecd <- out[, i, primfn] - maxPossiblePrim
      toSecd[toSecd < 0] <- 0
      getItems(toSecd, 3) <- secdfn
      out[, i, secdfn] <- out[, i, secdfn] + toSecd
      out[, i, primfn] <- mpmin(out[, i, primfn], maxPossiblePrim)

      woodland <- out[, , getItems(harvestAgg, 3)]
      stopifnot(harvestAgg <= woodland[, i - 1, ] / timestepLength,
                woodland[, i, primfn] <= woodland[, i - 1, primfn] - timestepLength * harvestAgg[, , primfn])
    }

    # consistency checks wood harvest area
    toolExpectLessDiff(harvest[, histYears, ], harvestHist, 0,
                       "In historical period, wood harvest area was not changed")
    toolExpectTrue(min(harvest) >= 0, "wood harvest area is >= 0")

    histYears <- getYears(harvest, as.integer = TRUE)
    histYears <- histYears[histYears < transitionYears[1]]
    toolCheckWoodHarvestArea(harvest[, histYears, ], out[, histYears, ],
                             "In historical period, ")

    futureYears <- setdiff(getYears(harvest, as.integer = TRUE), histYears)
    toolCheckWoodHarvestArea(harvest[, futureYears, ], out[, futureYears, ],
                             "After historical period, ")
  }

  # consistency checks land
  toolExpectLessDiff(out[, histYears, ], xTarget, 0,
                     "In historical period, land was not changed")
  toolExpectLessDiff(dimSums(out, 3), targetArea, 10^-5, "Total area is constant over time")
  toolPrimExpansionCheck(out)

  return(list(x = out,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Extrapolated land target data for harmonization",
              woodHarvestArea = harvest))
}
