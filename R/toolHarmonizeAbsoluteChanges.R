#' toolHarmonizeAbsoluteChanges
#'
#' Harmonize datasets by applying the absolute changes of the input data
#' relative to a single harmonization year on top of the target data at that
#' year. This is the aneris/mip "offset" method without convergence, i.e. a
#' permanent offset that is never faded out.
#'
#' @param xInput magpie object with the data to take absolute changes from,
#' usually model projections
#' @param xTarget magpie object with the historical target data
#' @param harmonizationPeriod a single integer year that must exist in both
#' xInput and xTarget; before this year the target dataset is used, from this
#' year onward the target value at harmonizationPeriod plus the absolute
#' change of the input data (relative to harmonizationPeriod) is used
#' @param level passed to toolReplaceExpansion and (offset by one) to the
#' toolStatusMessage calls in this function, see \code{\link{toolStatusMessage}}
#' @return magpie object with harmonized data. If a forest category (pltns,
#' primf, secdf) would go negative, the shortfall is instead funded from
#' other land (primn, secdn), deducted proportionally to their current
#' shares and added back to the forest categories proportionally to their
#' shortfall; a warning is emitted reporting the number and percentage of
#' affected cells and the min, mean and median of the negative values. Any
#' remaining negative values (forest cells where primn/secdn were
#' insufficient, and all non-forest negatives) are clamped to 0, again with
#' a warning if any clamping happened. The resulting surplus is then
#' redistributed: it is taken proportionally from the categories that still
#' have a positive value, excluding primf and primn, so the absolute-change
#' signal of primf/primn is preserved exactly. If the surplus cannot be
#' fully absorbed this way (e.g. because the eligible categories do not have
#' enough headroom), the previous behavior of renormalizing all categories
#' of the affected cell/timestep proportionally is used as a fallback,
#' emitting a note reporting how often this happened. Each cell/timestep
#' total area is kept constant throughout.
#' @author Pascal Sauer
toolHarmonizeAbsoluteChanges <- function(xInput, xTarget, harmonizationPeriod, level = 3) {
  y <- harmonizationPeriod

  stopifnot(!anyNA(xInput), !anyNA(xTarget))

  inputYears <- getYears(xInput, as.integer = TRUE)
  targetYears <- getYears(xTarget, as.integer = TRUE)

  stopifnot(length(y) == 1,
            round(y) == y,
            y %in% targetYears,
            y %in% inputYears,
            setequal(getItems(xInput, dim = 1), getItems(xTarget, dim = 1)),
            setequal(getItems(xInput, dim = 3), getItems(xTarget, dim = 3)))
  xInput <- xInput[getItems(xTarget, 1), , getItems(xTarget, 3)]

  futureYears <- inputYears[inputYears > y]
  stopifnot(length(futureYears) >= 1)
  delta <- xInput[, futureYears, ] - setYears(xInput[, y, ], NULL)
  changed <- setYears(xTarget[, y, ], NULL) + delta

  # forest rescue: instead of letting forest categories go negative, fund the shortfall
  # from other land (primn/secdn), deducted proportionally to their current shares and
  # added back to the forest categories proportionally to their shortfall. Remaining
  # (unfunded) negatives, as well as all non-forest negatives, are handled by the
  # clamping step below. Categories can legitimately be absent depending on the
  # recategorization, hence the intersect() guards.
  forest <- intersect(c("pltns", "primf", "secdf"), getItems(changed, dim = 3))
  otherLand <- intersect(c("primn", "secdn"), getItems(changed, dim = 3))
  if (length(forest) > 0 && any(changed[, , forest] < 0)) {
    forestValues <- as.vector(changed[, , forest])
    negativeForestValues <- forestValues[forestValues < 0]
    nNegativeForest <- length(negativeForestValues)
    shareNegativeForest <- nNegativeForest / length(forestValues)
    toolStatusMessage("warn",
                      paste0("funding ", nNegativeForest, " negative forest values (",
                             signif(100 * shareNegativeForest, 3),
                             "% of the forest values after the harmonization year) from other land, ",
                             "min value: ", signif(min(negativeForestValues), 3),
                             ", mean: ", signif(mean(negativeForestValues), 3),
                             ", median: ", signif(median(negativeForestValues), 3)),
                      level = level - 1)

    if (length(otherLand) > 0) {
      # shortfall per forest category/cell/timestep, and the amount that can actually
      # be funded per cell/timestep (limited by how much other land is available)
      shortfall <- pmax(-changed[, , forest], 0)
      shortfallTotal <- dimSums(shortfall, dim = 3)
      otherLandTotal <- pmax(dimSums(changed[, , otherLand], dim = 3), 0)
      funded <- pmin(shortfallTotal, otherLandTotal)

      # deduct the funded amount from primn/secdn proportionally to their current shares
      otherShare <- changed[, , otherLand] / otherLandTotal
      otherShare[!is.finite(otherShare)] <- 0 # cell/timestep without any other land
      changed[, , otherLand] <- changed[, , otherLand] - otherShare * funded

      # add the funded amount back to the forest categories that were negative,
      # proportionally to their shortfall (relevant if more than one is negative)
      forestShare <- shortfall / shortfallTotal
      forestShare[!is.finite(forestShare)] <- 0 # cell/timestep without any shortfall
      changed[, , forest] <- changed[, , forest] + forestShare * funded
    }
  }

  # clamping/redistribution only ever touches "changed" (the years after the
  # harmonization year); this keeps the years <= y bit-for-bit identical to
  # xTarget, which is not just cleaner but also avoids reintroducing floating
  # point noise into years that never needed any correction
  out <- mbind(xTarget[, targetYears <= y, ], changed)

  # only "changed" can contain negative values, xTarget is checked to be >= 0 where it is created
  if (any(changed < 0)) {
    nNegative <- sum(changed < 0)
    shareNegative <- nNegative / length(changed)
    minValue <- min(changed)
    toolStatusMessage("warn",
                      paste0("clamping ", nNegative, " negative values to 0 (",
                             signif(100 * shareNegative, 3), "% of the values after the harmonization year",
                             ", min value: ", signif(minValue, 3), ")"),
                      level = level - 1)
    changed[changed < 0] <- 0

    # redistribute the surplus created by clamping among the categories that
    # still have headroom, excluding primf/primn, so their absolute-change
    # signal is preserved exactly
    targetArea <- setYears(dimSums(xTarget[, y, ], dim = 3), NULL)
    # surplus is conceptually always >= 0 (clamping can only increase the cell/timestep
    # total), clip away tiny negative floating point noise to avoid spurious
    # "eligibleSum > surplus" matches (which would divide by zero if eligibleSum == 0)
    surplus <- pmax(dimSums(changed, dim = 3) - targetArea, 0)

    eligible <- setdiff(getItems(changed, dim = 3), c("primf", "primn"))
    eligibleSum <- dimSums(changed[, , eligible], dim = 3)

    scaling <- eligibleSum
    scaling[] <- 1
    absorbable <- eligibleSum > surplus
    scaling[absorbable] <- (eligibleSum[absorbable] - surplus[absorbable]) / eligibleSum[absorbable]
    changed[, , eligible] <- changed[, , eligible] * scaling

    # fallback for cell/timesteps where the eligible categories could not
    # absorb the full surplus (surplus >= eligibleSum, including
    # eligibleSum == 0); this reproduces the old proportional renormalization,
    # but only where it is actually needed, so everywhere else the
    # absolute-change signal is returned without any renormalization
    needsFallback <- !absorbable & surplus > 0
    if (any(needsFallback)) {
      normalization <- targetArea / dimSums(changed, dim = 3)
      normalization[!needsFallback] <- 1
      normalization[!is.finite(normalization)] <- 1 # cell/timestep without any land
      toolStatusMessage("note",
                        paste0(sum(needsFallback), " cell/timesteps needed the fallback renormalization ",
                               "(surplus could not be fully absorbed by the eligible categories)"),
                        level = level - 1)
      changed <- changed * normalization
    }

    out <- mbind(xTarget[, targetYears <= y, ], changed)
  }

  out <- toolReplaceExpansion(out, "primf", "secdf", warnThreshold = 100, level = level)
  out <- toolReplaceExpansion(out, "primn", "secdn", warnThreshold = 100, level = level)

  return(out)
}
