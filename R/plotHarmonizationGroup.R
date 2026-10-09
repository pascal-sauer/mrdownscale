#' plotHarmonizationGroup
#'
#' Internal helper (not exported) that plots one harmonization group: global
#' totals per category and harmonization method as a line chart, written as a
#' PNG to disk.
#'
#' @param prefix filename prefix distinguishing the data type, e.g. "land" or "nonland"
#' @param name group name, used for the plot title and the output filename
#' @param d data.frame with columns Year, Item, Value and Method, holding mapped
#' items that are not yet aggregated (they are summed per Year/Item/Method here)
#' @param ylab y axis label, e.g. "area [Mha]"
#' @param harmonizationPeriod two integer years, marked by vertical lines (the
#' second one only if a fade harmonization is present)
#' @param plotdir directory to write the PNG to, by default the working directory
#' @return Invisibly, the path of the written PNG file.
#' @author Pascal Sauer
#' @seealso \code{\link{plotHarmonizedLand}}, \code{\link{plotHarmonizedNonLand}}
plotHarmonizationGroup <- function(prefix, name, d, ylab, harmonizationPeriod, plotdir = ".") {
  linetypes <- c(absoluteChanges = "solid", inputUnharmonized = "dashed",
                 fadeForest = "dotted", fade = "dotdash")
  fontsize <- 14
  d <- aggregate(Value ~ Year + Item + Method, d, sum)
  d$Item <- factor(d$Item)
  fadeEndShown <- any(c("fade", "fadeForest") %in% unique(d$Method))
  .data <- ggplot2::.data
  p <- ggplot2::ggplot(d, ggplot2::aes(.data$Year, .data$Value, color = .data$Item, linetype = .data$Method)) +
    ggplot2::theme_bw(base_size = fontsize) +
    ggplot2::geom_line() +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(harmonizationPeriod[1],
                                       if (fadeEndShown) harmonizationPeriod[2]), linewidth = 0.3) +
    ggplot2::scale_linetype_manual(values = linetypes[unique(d$Method)],
                                   name = NULL,
                                   guide = ggplot2::guide_legend(keywidth = 3)) +
    ggplot2::labs(title = paste0("harmonization comparison, group ", name),
                  y = ylab, color = "category") +
    ggplot2::theme(legend.text = ggplot2::element_text(size = fontsize - 3))
  file <- file.path(plotdir, paste0("plot-harmonized-", prefix, "-", name, ".png"))
  ggplot2::ggsave(file, p, width = 10, height = 6, dpi = 120)
  message("wrote ", file)
  return(invisible(file))
}
