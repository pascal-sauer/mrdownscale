#' toolReplaceExpansion
#'
#' Expansion from one timestep to the next of one land type is replaced with
#' another land type.
#'
#' @param x a magpie object
#' @param from name of a land category, e.g. "primf"
#' Expansion of this category happening from one timestep to the next
#' will be replaced. If missing in x return x unchanged
#' @param to name of another land category, e.g. "secdf"
#' Expansion of 'from' will be replaced with expansion of 'to'
#' @param ... not used, will throw an error if supplied
#' @param noteThreshold expansion greater than this will trigger a note
#' that expansion was replaced
#' @param warnThreshold expansion greater than this will trigger a warning
#' that expanding considerably
#' @return a magpie object with expansion of 'from' replaced by 'to'
#' @author Pascal Sauer
#' @export
toolReplaceExpansion <- function(x, from, to, ..., noteThreshold = 10^-10, warnThreshold = 10^-5) {
  if (!from %in% getItems(x, dim = 3)) {
    return(x)
  }
  stopifnot(nyears(x) >= 2,
            from != to,
            to %in% getItems(x, dim = 3),
            ...length() == 0,
            noteThreshold <= warnThreshold)

  maxDiff <- 0
  for (i in 2:nyears(x)) {
    difference <- x[, i, from] - setYears(x[, i - 1, from], getYears(x)[i])
    difference[difference < 0] <- 0
    maxDiff <- max(maxDiff, max(difference))
    # use pmin instead of subtracting to avoid tiny expansions due to numerical imprecision
    x[, i, from] <- pmin(x[, i, from], setYears(x[, i - 1, from], getYears(x)[i]))
    x[, i, to] <- x[, i, to] + magclass::setNames(difference, to)
  }

  if (maxDiff > warnThreshold) {
    toolStatusMessage("warn", paste0(from, " is expanding considerably"))
  }

  if (maxDiff > noteThreshold) {
    toolStatusMessage("note", paste0("replaced ", from, " expansion ",
                                     "(max expansion: ", signif(maxDiff, 3), ") with ", to, " expansion"))
  }

  return(x)
}
