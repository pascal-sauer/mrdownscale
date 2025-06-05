#' toolGetHarmonizer
#'
#' Get a harmonizer function by name.
#'
#' @param harmonizerName name of a harmonizer function, currently offset, fade, fadeForest
#' @return harmonizer function
#' @seealso \code{\link{toolHarmonizeOffset}}, \code{\link{toolHarmonizeFade}}, \code{\link{toolHarmonizeFadeForest}}
#' @author Pascal Sauer
toolGetHarmonizer <- function(harmonizerName) {
  # function(...) toolHarmonizeOffset(...) instead of passing
  # toolHarmonizeOffset directly so madrat recognizes it as dependency
  harmonizers <- list(offset = function(...) toolHarmonizeOffset(...),
                      fade = function(...) toolHarmonizeFade(...),
                      fadeForest = function(...) toolHarmonizeFadeForest(...))
  stopifnot(harmonizerName %in% names(harmonizers))
  return(harmonizers[[harmonizerName]])
}
