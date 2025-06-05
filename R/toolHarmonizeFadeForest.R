#' toolHarmonizeFadeForest
#'
#' Harmonize using \code{\link{toolHarmonizeFade}}, then sum up primf + secdf,
#' then disaggregate again in a way that minimizes primf to secdf conversion.
#'
#' @inheritParams toolHarmonizeFade
#' @inherit toolHarmonizeFade return
#' @author Pascal Sauer
toolHarmonizeFadeForest <- function(xInput, xTarget, harmonizationPeriod) {
  x <- toolHarmonizeFade(xInput, xTarget, harmonizationPeriod, level = 4)
  years <- getYears(x, as.integer = TRUE)
  psf <- c("primf", "secdf")
  if (!all(psf %in% getItems(x))) {
    return(x)
  }
  hp1 <- harmonizationPeriod[1]
  forest <- dimSums(x[, years >= hp1, psf], 3)

  primf <- xTarget[, hp1, "primf"]
  for (year in years[years > hp1]) {
    previousYear <- primf[, nyears(primf), ]
    p <- mpmin(previousYear, forest[, year, ])
    p <- setYears(p, year)
    primf <- mbind(primf, p)
  }
  primsecdf <- mbind(primf,
                     magclass::setNames(forest - primf, "secdf"))

  stopifnot(primsecdf >= 0,
            abs(forest - dimSums(primsecdf, 3)) < 10^-10)
  x[, years >= hp1, psf] <- primsecdf
  return(x)
}
