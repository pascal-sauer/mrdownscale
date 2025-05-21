toolEqualizeArea <- function(xInput, xTarget) {
  inSum <- dimSums(xInput, dim = 3)
  tSum <- dimSums(xTarget, dim = 3)
  if (max(abs(inSum[, 1, ] - tSum[, 1, ])) >= 10^-5) {
    corr <- setYears(dimSums(xTarget[, 1, ], dim = 3) / dimSums(xInput[, 1, ], dim = 3), NULL)
    stopifnot(is.finite(corr), corr >= 0)
    toolStatusMessage("note", paste0("input data multiplied with correction factors to match target areas ",
                                     "(max ratio = ", round(max(corr), 2),
                                     ", min ratio = ", round(min(corr), 2),  ")"))
    return(xInput * corr)
  }
  return(xInput)
}
