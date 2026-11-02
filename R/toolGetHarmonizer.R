#' toolGetHarmonizer
#'
#' Get a harmonizer function by name.
#'
#' @param harmonizerName name of a harmonizer function, currently offset, fade, fadeForest,
#' absoluteChanges
#' @param constantTotal TRUE for area conservative land data, FALSE otherwise, only
#' used by the absoluteChanges harmonizer, ignored by the others
#' @return harmonizer function
#' @seealso \code{\link{toolHarmonizeOffset}}, \code{\link{toolHarmonizeFade}},
#' \code{\link{toolHarmonizeFadeForest}}, \code{\link{toolHarmonizeAbsoluteChanges}}
#' @author Pascal Sauer
toolGetHarmonizer <- function(harmonizerName, constantTotal) {
  # function(...) toolHarmonizeOffset(...) instead of passing
  # toolHarmonizeOffset directly so madrat recognizes it as dependency
  harmonizers <- list(offset = function(...) toolHarmonizeOffset(...),
                      fade = function(...) toolHarmonizeFade(...),
                      fadeForest = function(...) toolHarmonizeFadeForest(...),
                      absoluteChanges = function(...) toolHarmonizeAbsoluteChanges(..., constantTotal = constantTotal))
  stopifnot(harmonizerName %in% names(harmonizers))
  return(harmonizers[[harmonizerName]])
}
