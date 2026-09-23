#' toolPrimExpansionCheck
#'
#' Check that primary forest (primf) and primary nature/nonland (primf) are
#' never expanding using toolExpectTrue. If they are expanding, a warning
#' including the maximum expansion is thrown. Categories not present
#' in x are skipped.
#'
#' @param x A magclass object
#'
#' @author Pascal Sauer
toolPrimExpansionCheck <- function(x) {
  stopifnot(!anyNA(x))
  for (v in intersect(getItems(x, 3), c("primf", "primn"))) {
    maxExpansion <- toolMaxExpansion(x[, , v])
    toolExpectTrue(maxExpansion <= 0,
                   paste0(v, " is never expanding",
                          if (maxExpansion > 0) paste0(", max expansion: ", signif(maxExpansion, 3))),
                   falseStatus = "warn")
  }
}
