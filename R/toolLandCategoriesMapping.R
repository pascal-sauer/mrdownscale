#' toolLandCategoriesMapping
#'
#' Computes a land category mapping between input and target data set via
#' retrieving mappings from these sources to an internally defined
#' reference categorization and merging them to a direct input-to-target
#' mapping.
#' Mappings to reference categories are stored internally in the package and
#' have to be added to it if new input and/or target data should be supported.
#'
#' @param input name of the land input source to be used
#' @param target name of the land target source to be used
#' @author Jan Philipp Dietrich
#'
toolLandCategoriesMapping <- function(input, target) {
  .getMap <- function(x) {
    if (x == "magpie") {
      out <- toolGetMapping("magpie2ref.csv", where = "mrdownscale")
    } else if (x == "landuseinit") {
      out <- toolGetMapping("landuseinit2ref.csv", where = "mrdownscale")
    } else if (x == "luh2") {
      out <- toolGetMapping("luh2ref.csv", where = "mrdownscale")
    } else if (x == "luh2mod") {
      out <- toolGetMapping("luhmod2ref.csv", where = "mrdownscale")
    } else {
      stop("Categories mapping for type \"", x, "\" not available!", call. = FALSE)
    }
    return(out)
  }
  input2ref  <- .getMap(input)
  output2ref <- .getMap(target)

  if (!setequal(input2ref$reference, output2ref$reference)) {
    warning("Input map and output map contain inconsistent reference information")
  }
  map <- merge(input2ref, output2ref, by = "reference", suffixes = c("Input", "Output"))
  map$merge <- paste0(map$dataInput, "__", map$dataOutput)
  if (anyDuplicated(map$reference)) {
    warning("Insuficient granularity of reference categories, as a reference category is mapped more than once (\"",
            paste(unique(map$reference[duplicated(map$reference)]), collapse = "\", \""), "\").")
  }
  return(map)
}
