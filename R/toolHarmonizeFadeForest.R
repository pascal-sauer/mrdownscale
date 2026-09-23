#' toolHarmonizeFadeForest
#'
#' Harmonize using \code{\link{toolHarmonizeFade}}, then sum up primf + secdf,
#' then disaggregate again in a way that minimizes primf to secdf conversion,
#' without exceeding linear extrapolation of primf.
#'
#' @inheritParams toolHarmonizeFade
#' @inherit toolHarmonizeFade return
#' @author Pascal Sauer
toolHarmonizeFadeForest <- function(xInput, xTarget, harmonizationPeriod) {
  hp1 <- harmonizationPeriod[1]

  x <- toolHarmonizeFade(xInput, xTarget, harmonizationPeriod)
  years <- getYears(x, as.integer = TRUE)
  psf <- c("primf", "secdf")
  if (!all(psf %in% getItems(x, 3))) {
    return(x)
  }
  forest <- dimSums(x[, years >= hp1, psf], 3)

  # linear extrapolation of primf, don't want to exceed this, also ensures some primf is always harvested
  primfTarget <- xTarget[, getYears(xTarget, TRUE) <= hp1, "primf"]
  decline <- -primfTarget[, -1, ] + setYears(primfTarget[, -nyears(primfTarget), ],
                                             getYears(primfTarget[, -1, ]))
  primfMeanDecline <- magpply(decline, mean, DIM = 2)
  linearEx <- xTarget[, hp1, "primf"]
  for (year in years[years > hp1]) {
    linearEx <- mbind(linearEx, setYears(linearEx[, nyears(linearEx), ] - primfMeanDecline, year))
  }
  linearEx[linearEx < 0] <- 0

  primf <- xTarget[, hp1, "primf"]
  for (year in years[years > hp1]) {
    p <- linearEx[, year, ]
    p <- pmin(p,
              magclass::setNames(
                forest[, year, ],
                "primf"
              ),
              setYears(
                primf[, nyears(primf), ],
                year
              ))
    primf <- mbind(primf, p)
  }
  primsecdf <- mbind(primf,
                     magclass::setNames(forest - primf, "secdf"))

  stopifnot(primsecdf >= 0,
            abs(forest - dimSums(primsecdf, 3)) < 10^-10)
  x[, years >= hp1, psf] <- primsecdf
  return(x)
}
