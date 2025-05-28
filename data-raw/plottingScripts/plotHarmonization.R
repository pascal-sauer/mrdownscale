# nolint start
plotHarmonization <- function(input, target, harmonizationPeriod) {
  lir <- calcOutput("LandInputRecategorized", input = input, target = target, aggregate = FALSE)
  ltex <- calcOutput("LandTargetExtrapolated", input = input, target = target,
                     harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
  lh <- calcOutput("LandHarmonized", aggregate = FALSE, input = input,
                   target = target, harmonizationPeriod = harmonizationPeriod)

  crops <- grep("^c[34](ann|per|nfx)", getItems(lir, 3), value = TRUE)
  years <- getYears(lh, TRUE)
  years <- years[years >= harmonizationPeriod[[1]] & years < harmonizationPeriod[[2]]]

  lir <- dimSums(lir[, years, crops], c(1, 3))
  ltex <- dimSums(ltex[, years, crops], c(1, 3))
  lh <- dimSums(lh[, years, crops], c(1, 3))
  x <- mbind(magclass::setNames(lir, "model"),
            #  magclass::setNames(ltex, "historical extrapolated"),
             magclass::setNames(lh, "harmonized"))
  x <- as.data.frame(x)
  names(x)[4:5] <- c("Data", "CroplandMha")
  withr::local_package("ggplot2")
  fontsize <- 25
  ggplot(data = x, aes(x = Year, y = CroplandMha, group = Data, color = Data)) +
    geom_line(size = 1.5) +
    theme(text = element_text(size = fontsize),
          axis.text = element_text(size = fontsize),
          axis.title = element_text(size = fontsize),
          plot.title = element_text(size = fontsize),
          legend.text = element_text(size = fontsize),
          legend.title = element_text(size = fontsize)) +
    scale_y_continuous(n.breaks = 8)
}

plotHistogram <- function(input, target, harmonizationPeriod, variables) {
  # e.g. crops: variables <- grep("^c[34](ann|per|nfx)", getItems(lir, 3), value = TRUE)
  lir <- calcOutput("LandInputRecategorized", input = input, target = target, aggregate = FALSE)
  ltex <- calcOutput("LandTargetExtrapolated", input = input, target = target,
                     harmonizationPeriod = harmonizationPeriod, aggregate = FALSE)
  lh <- calcOutput("LandHarmonized", aggregate = FALSE, input = input,
                   target = target, harmonizationPeriod = harmonizationPeriod)

  years <- getYears(lh, TRUE)
  years <- years[years >= harmonizationPeriod[[1]] & years < harmonizationPeriod[[2]]]

  lirCrop <- dimSums(lir[, years, variables], 3)
  ltexCrop <- dimSums(ltex[, years, variables], 3)
  lhCrop <- dimSums(lh[, years, variables], 3)

  x <- lirCrop - lhCrop
  x <- x / dimSums(lir[, years, ], 3)

  x <- as.data.frame(x)
  # x <- x[x$Region == "EUR", ]
  names(x)[4:5] <- c("Data", "RelativeDiff")
  withr::local_package("ggplot2")
  fontsize <- 25
  (ggplot(data = x, aes(RelativeDiff, group = Region))
    + geom_histogram()
    + facet_grid(vars(Year), vars(Region))
    # + facet_wrap(~Year)
    + theme(text = element_text(size = fontsize),
            axis.text = element_text(size = fontsize),
            axis.title = element_text(size = fontsize),
            plot.title = element_text(size = fontsize),
            legend.text = element_text(size = fontsize),
            legend.title = element_text(size = fontsize))
    # + scale_y_continuous(n.breaks = 5)
    # + scale_x_continuous(n.breaks = 8)
  )
}
# nolint end
