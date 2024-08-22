#' calcLandTargetExtrapolated
#'
#' First aggregate to low resolution, then extrapolate to the given years
#' using toolExtrapolate.
#'
#' @param target name of the target dataset, options are: luh2, luh2mod
#' luh2mod will split secdf into forestry and secdf
#' @return extrapolated land target data
#' @author Pascal Sauer
calcLandTargetExtrapolated <- function(input = "magpie", target = "luh2mod",
                                       harmonizationPeriod = c(2015, 2050), extrapolate = TRUE) {
  xInput <- calcOutput("LandHarmonizedCategories", input = input,
                       target = target, aggregate = FALSE)
  xTarget <- calcOutput("LandTarget", target = target, aggregate = FALSE)

  # bring target data to spatial resolution of input data
  ref    <- as.SpatVector(xInput[, 1, 1])[, c(".region", ".id")]
  xTarget <- terra::extract(xTarget, ref, sum, na.rm = TRUE, bind = TRUE)
  xTarget <- as.magpie(xTarget)
  stopifnot(setequal(getItems(xInput, 3), getItems(xTarget, 3)))
  xTarget <- xTarget[, , getItems(xInput, 3)] # harmonize order of dim 3

  if (extrapolate) {
    # ---------- extrapolate -------------
    inputYears <- getYears(xInput, as.integer = TRUE)
    transitionYears <- inputYears[inputYears > harmonizationPeriod[1] & inputYears < harmonizationPeriod[2]]
    exTarget <- toolExtrapolate(xTarget, transitionYears)
    exTarget[exTarget < 0] <- 0

    # ---------- normalize -------------
    # normalize exTarget so that its total sum over all layers agrees for all time steps
    # with the sum over all layers in target in the harmonization year (e.g. makes sure
    # that the land changes in a land data set do not alter the total sum of land.)
    targetArea <- dimSums(setYears(xTarget[, harmonizationPeriod[1], ], NULL), dim = 3)
    exTarget <- exTarget * targetArea / dimSums(exTarget, dim = 3)
    exTarget[is.na(exTarget)] <- 0

    # ------- calculate wood harvest shares -------
    harvest <- calcOutput("NonlandLowRes", input = input, target = target, aggregate = FALSE)
    harvest <- harvest[, , endsWith(getItems(harvest, 3), "wood_harvest_area")]
    harvest <- .aggregateWoodHarvest(harvest)
    timestepLength <- unique(diff(getYears(harvest, as.integer = TRUE)))
    stopifnot(length(timestepLength) == 1)
    harvest <- harvest * timestepLength # harvest is per year, need per timestep
    stopifnot(identical(dimnames(harvest)[1:2], dimnames(xTarget)[1:2]))

    # calculate share: primf|primn wood harvest area / total primf|primn area
    primShare <- dimSums(harvest[, , c("primf", "primn")], 2) / dimSums(xTarget[, , c("primf", "primn")], 2)
    primShare[is.na(primShare)] <- 0
    primShare[primShare > 1] <- 1 # TODO note this in the log
    stopifnot(0 <= primShare, primShare <= 1)

    # calculate share: forest (primf + secdf) wood harvest area / total forest area
    forest <- c("primf", "secdf")
    totalShareForest <- dimSums(harvest[, , forest], c(2, 3)) / dimSums(xTarget[, , forest], c(2, 3))
    totalShareForest[is.na(totalShareForest)] <- 0
    totalShareForest[totalShareForest > 1] <- 1 # TODO note this in the log
    stopifnot(0 <= totalShareForest, totalShareForest <= 1)

    # calculate share: nature (primn + secdn) wood harvest area / total nature area
    nature <- c("primn", "secdn")
    totalShareNature <- dimSums(harvest[, , nature], c(2, 3)) / dimSums(xTarget[, , nature], c(2, 3))
    totalShareNature[is.na(totalShareNature)] <- 0
    totalShareNature[totalShareNature > 1] <- 1 # TODO note this in the log
    stopifnot(0 <= totalShareNature, totalShareNature <= 1)

    out <- mbind(xTarget, exTarget)

    outHarvest <- add_columns(harvest, paste0("y", transitionYears), dim = 2)

    for (i in match(transitionYears, getYears(out, as.integer = TRUE))) {
      year <- getYears(out, as.integer = TRUE)[i]
      maxPossiblePrimf <- out[, i - 1, "primf"] - outHarvest[, i - 1, "primf"]
      stopifnot(maxPossiblePrimf >= 0)
      toSecdf <- out[, i, "primf"] - maxPossiblePrimf
      toSecdf[toSecdf < 0] <- 0
      out[, i, "secdf"] <- out[, i, "secdf"] + toSecdf
      out[, i, "primf"] <- pmin(out[, i, "primf"], maxPossiblePrimf)

      outHarvest[, i, "primf"] <- outHarvest[, i, "primf"] * primShare[, , "primf"]
      totalHarvest <- dimSums(out[, i, c("primf", "secdf")], 3) * totalShareForest
      outHarvest[, i, "secdf"] <- totalHarvest - outHarvest[, i, "primf"]
      stopifnot(outHarvest[, i, "secdf"] >= 0)
      if (outHarvest[, i, "secdf"] > out[, i, "secdf"]) {
        browser() # TODO
      }
    }
    stopifnot(is.finite(outHarvest), outHarvest >= 0)
    # TODO same for primn
    outHarvest <- outHarvest / timestepLength
  } else {
    out <- xTarget
  }

  # TODO checks

  return(list(x = out,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Extrapolated land target data for harmonization"))
}
