#' toolReportAreaDeviation
#'
#' Reports, per category group and per variable, how much the corrected output
#' deviates from the raw input data: the positive and the negative deviations
#' (xOut - xRaw) accumulated separately over all cells and timesteps, as mean
#' per timestep in Mha and as a percentage of the output area accumulated over
#' the same timesteps, plus the net change as the sum of both columns, in six
#' columns: =%, +%, -%, =Mha, +Mha, -Mha.
#'
#' @param xRaw harmonized data before corrections, as magpie object, with the
#'   years of xOut after the harmonization year
#' @param xOut harmonized data after corrections, as magpie object, including
#'   the harmonization year and the years before it
#' @param groups Named list of character vectors with the category groups to be
#'   reported separately, named as they should appear in the report; every
#'   category must be in one of the groups or be "urban"
#'
#' @author Pascal Sauer
toolReportAreaDeviation <- function(xRaw, xOut, groups) {
  stopifnot(identical(getItems(xRaw, 3), getItems(xOut, 3)))
  yearsAfter <- getYears(xRaw, as.integer = TRUE)
  deviation <- xOut[, yearsAfter, ] - xRaw
  .accumulateByItem <- function(x) {
    result <- as.numeric(dimSums(x, c(1, 2)))
    result[abs(result) < 10^-5] <- 0
    names(result) <- getItems(xOut, 3)
    return(result)
  }
  positive <- .accumulateByItem(pmax(deviation, 0))
  negative <- .accumulateByItem(pmin(deviation, 0))
  net <- positive + negative
  outputArea <- .accumulateByItem(xOut[, yearsAfter, ])

  .formatNumber <- function(x) as.character(signif(x, 3))
  .meanPerTimestep <- function(x) x / length(yearsAfter)
  # share of the accumulated output area, as percent
  .formatPct <- function(part, total) {
    pct <- signif(100 * part / total, 3)
    return(ifelse(is.na(pct), "n/a", as.character(pct)))
  }
  # align a column of number-like strings vertically at the decimal separator
  .alignDecimal <- function(x) {
    intPart <- sub("\\..*$", "", x)
    fracPart <- ifelse(grepl(".", x, fixed = TRUE), sub(".*\\.", "", x), NA_character_)
    fracWidth <- suppressWarnings(max(nchar(fracPart), na.rm = TRUE)) # suppress because all NA would warn
    intText <- format(intPart, width = max(nchar(intPart)), justify = "right")
    if (is.finite(fracWidth) && fracWidth > 0) {
      fracText <- ifelse(is.na(fracPart), "", fracPart)
      return(paste0(intText,
                    ifelse(is.na(fracPart), " ", "."),
                    fracText,
                    strrep(" ", fracWidth - nchar(fracText))))
    }
    return(intText)
  }

  for (groupName in names(groups)) {
    groupItems <- groups[[groupName]]
    if (length(groupItems) == 0) {
      next
    }
    groupPositive <- positive[groupItems]
    groupNegative <- negative[groupItems]
    groupNet <- net[groupItems]
    groupOutputArea <- outputArea[groupItems]
    # sort by the gross deviation, so the variables that deviate most surface
    # first, ties in group order
    orderIndex <- order(groupPositive - groupNegative, decreasing = TRUE)
    groupPositive <- groupPositive[orderIndex]
    groupNegative <- groupNegative[orderIndex]
    groupNet <- groupNet[orderIndex]
    groupOutputArea <- groupOutputArea[orderIndex]
    groupTotalPositive <- sum(groupPositive)
    groupTotalNegative <- sum(groupNegative)
    groupTotalNet <- groupTotalPositive + groupTotalNegative
    groupTotalOutputArea <- sum(groupOutputArea)
    tableRows <- rbind(
      c("variable", "=%", "+%", "-%", "=Mha", "+Mha", "-Mha"),
      c("total", .formatPct(groupTotalNet, groupTotalOutputArea),
        .formatPct(groupTotalPositive, groupTotalOutputArea),
        .formatPct(groupTotalNegative, groupTotalOutputArea),
        .formatNumber(.meanPerTimestep(groupTotalNet)),
        .formatNumber(.meanPerTimestep(groupTotalPositive)),
        .formatNumber(.meanPerTimestep(groupTotalNegative))),
      cbind(names(groupPositive),
            .formatPct(groupNet, groupOutputArea),
            .formatPct(groupPositive, groupOutputArea),
            .formatPct(groupNegative, groupOutputArea),
            .formatNumber(.meanPerTimestep(groupNet)),
            .formatNumber(.meanPerTimestep(groupPositive)),
            .formatNumber(.meanPerTimestep(groupNegative)))
    )
    tableText <- apply(cbind(format(tableRows[, 1], width = max(nchar(tableRows[, 1]))),
                             .alignDecimal(tableRows[, 2]),
                             .alignDecimal(tableRows[, 3]),
                             .alignDecimal(tableRows[, 4]),
                             .alignDecimal(tableRows[, 5]),
                             .alignDecimal(tableRows[, 6]),
                             .alignDecimal(tableRows[, 7])),
                       1, paste, collapse = ", ")
    toolStatusMessage("note",
                      paste0("Deviation output vs. history plus raw input trend in group \"", groupName,
                             "\" (deviations accumulated over cells and timesteps, as mean per timestep",
                             " in Mha and as % of the output area of the same period):\n",
                             paste(tableText, collapse = "\n")))
  }
  return(invisible(NULL))
}
