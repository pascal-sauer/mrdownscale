toolEqualizeArea <- function(x, y, level = 1) {
  stopifnot(nyears(y) == 1)
  xSum <- dimSums(x[, 1, ], dim = 3)
  ySum <- dimSums(y, dim = 3)
  if (max(abs(xSum - ySum)) >= 10^-5) {
    corr <- setYears(ySum / xSum, NULL)
    corr[ySum == 0 & xSum == 0] <- 1
    stopifnot(is.finite(corr), corr >= 0)
    toolStatusMessage("note", paste0("correction factors were applied to equalize area ",
                                     "(max ratio = ", round(max(corr), 2),
                                     ", min ratio = ", round(min(corr), 2),  ")"),
                      level = level)
    x <- x * corr
    xSum <- dimSums(x[, 1, ], dim = 3)
  }

  stopifnot(max(abs(xSum - ySum)) < 10^-5)
  return(x)
}
