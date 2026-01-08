#' toolLandCategoriesMapping
#'
#' Computes a land category mapping between input and target data set via
#' retrieving mappings from these sources to an internally defined
#' reference categorization and merging them to a direct input-to-target
#' mapping.
#' Mappings to reference categories are stored internally in the package and
#' have to be added to it if new input and/or target data should be supported.
#'
#' @inheritParams calcLandInput
#' @param target name of the land target source to be used
#' @author Jan Philipp Dietrich
#'
toolLandCategoriesMapping <- function(input, target) {
  .getMap <- function(x) {
    if (x == "magpie") {
      out <- toolGetMapping("referenceMappings/magpie.csv", where = "mrdownscale")
    } else if (x == "witch") {
      out <- toolGetMapping("referenceMappings/witch.csv", where = "mrdownscale")
    } else if (x == "coffee") {
      out <- toolGetMapping("referenceMappings/coffee.csv", where = "mrdownscale")
    } else if (x == "landuseinit") {
      out <- toolGetMapping("referenceMappings/landuseinit.csv", where = "mrdownscale")
    } else if (x == "landuseinitchina") {
      out <- toolGetMapping("referenceMappings/landuseinitchina.csv", where = "mrdownscale")
    } else if (x == "luh2") {
      out <- toolGetMapping("referenceMappings/luh2.csv", where = "mrdownscale")
    } else if (x == "luh3pltns") {
      out <- toolGetMapping("referenceMappings/luh3pltns.csv", where = "mrdownscale")
    } else if (x == "luh2mod") {
      out <- toolGetMapping("referenceMappings/luh2mod.csv", where = "mrdownscale")
    } else {
      stop("Categories mapping for type \"", x, "\" not available!", call. = FALSE)
    }
    return(out)
  }
  input2ref  <- .getMap(input)
  target2ref <- .getMap(target)

  if (!setequal(input2ref$reference, target2ref$reference)) {
    warning("Input map and output map contain inconsistent reference information")
  }
  map <- merge(input2ref, target2ref, by = "reference", suffixes = c("Input", "Output"))
  map$merge <- paste0(map$dataInput, "__", map$dataOutput)
  if (anyDuplicated(map$reference)) {
    warning("Insuficient granularity of reference categories, as a reference category is mapped more than once (\"",
            paste(unique(map$reference[duplicated(map$reference)]), collapse = "\", \""), "\").")
  }
  return(map)
}
