#' toolRedistributeBioh
#'
#' Move the bioh of clusters without wood harvest area to the clusters that still
#' have some. Their recalculated bioh (harvest area times harvest rate) is
#' zero, so rates times land cannot equal target bioh and would be stranded.
#' Per subcategory and year, the stranded bioh is added to clusters with
#' recalculated bioh > 0 in the same region proportional to their recalculated
#' bioh, falling back to the global pool if a region has no such clusters.
#' Bioh that cannot be moved anywhere is reported and dropped, so the
#' returned bioh is zero everywhere where the recalculated bioh is zero.
#'
#' @param bioh bioh in kg C yr-1
#' @param biohCalculated bioh in in kg C yr-1 recalculated from harvest rate
#'   times harvest area
#' @return bioh in kg C yr-1 with the stranded bioh moved to the receivers,
#'   zero on clusters without harvest area, totals preserved unless bioh had
#'   to be dropped
#'
#' @author Pascal Sauer
toolRedistributeBioh <- function(bioh, biohCalculated) {
  receiver <- dimSums(biohCalculated, 3) > 0
  stranded <- bioh - bioh * receiver
  out <- bioh * receiver

  regionStranded <- dimSums(stranded, 1.2)
  regionWeight <- dimSums(biohCalculated, 1.2)
  regionReceiver <- regionWeight > 0
  # only needed where regions have neither stranded bioh nor receivers,
  # which yields a share of zero either way
  regionWeight[regionWeight == 0] <- 1
  out <- out + regionStranded * regionReceiver / regionWeight * biohCalculated

  # regions without receivers hand their stranded bioh to the global pool
  globalWeight <- dimSums(biohCalculated, 1)
  globalReceiver <- globalWeight > 0
  globalStranded <- dimSums(regionStranded * (1 - regionReceiver), 1) * globalReceiver
  globalWeight[globalWeight == 0] <- 1
  out <- out + globalStranded / globalWeight * biohCalculated

  # bioh of a subcategory and year without receivers anywhere cannot be moved
  # without losing totals, it is reported and dropped
  unmovable <- stranded * (1 - regionReceiver) * (1 - globalReceiver)
  totalBioh <- as.numeric(dimSums(bioh, c(1, 2, 3)))
  dropped <- as.numeric(dimSums(unmovable, c(1, 2, 3)))
  moved <- as.numeric(dimSums(stranded, c(1, 2, 3))) - dropped
  affectedClusterYears <- sum(as.numeric(dimSums(stranded, 3) > 0))
  totalClusterYears <- length(getItems(bioh, 1)) * length(getYears(bioh))
  .formatPct <- function(part, total) paste0(signif(100 * part / total, 3), "%")
  if (affectedClusterYears == 0) {
    toolStatusMessage("ok", "no stranded bioh")
  } else {
    clusterYearsText <- paste0(affectedClusterYears, " of ", totalClusterYears, " cluster years (",
                               .formatPct(affectedClusterYears, totalClusterYears), ")")
    if (moved > 0) {
      toolStatusMessage("note",
                        paste0("redistributed stranded bioh of clusters without harvest area: ", signif(moved, 3),
                               " kg C yr-1 (", .formatPct(moved, totalBioh),
                               ", accumulated over clusters and years) from ", clusterYearsText))
    }
    if (dropped > 0) {
      toolStatusMessage("warn",
                        paste0("dropped stranded bioh without receivers: ", signif(dropped, 3), " kg C yr-1 (",
                               .formatPct(dropped, totalBioh), ", accumulated over clusters and years) in ",
                               clusterYearsText))
    }
  }

  return(out)
}
