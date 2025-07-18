#' calcLandCategorizationWeight
#'
#' Computes weights for a weighted category mapping and aggregates the weights
#' to the spatial resolution described by the arguments \code{geometry} and \code{crs}.
#'
#' @note This calc-function has a rather unusual shape in that sense that the arguments to
#' be provided are not simple configuration settings but rather relative complex. For this
#' kind of implementation it is typically advised to use tool functions. The reason that a
#' calc function is used in this particular case instead is, that the inputs of this
#' function usually do not change over a long period of time so that caching of results
#' becomes key for the overall performance of the data processing, which is available
#' for calc- but not tool-functions. If tool-functions might support caching in the future
#' as well as a conversion to a tool function might become a sensible option.
#'
#' @param map a map in form of a 2 column data.frame containing a mapping between
#' input categories and output categories for the data dimension of a magpie object
#' @param geometry the geometry of the magpie object for which the categories should
#' be mapped as given in the geometry attribute of a magpie object with
#' geometry information \code{attr(x, "geometry")}.
#' @param crs the coordinate reference system as returned by \code{attr(x, "crs")} from
#' a magpie object with coordinates information.
#' @author Jan Philipp Dietrich
calcLandCategorizationWeight <- function(map, geometry, crs) {
  .getTarget <- function(geometry, crs) {
    target <- new.magpie(names(geometry), sets = c("id", "temporal", "data"))
    attr(target, "geometry") <- geometry
    attr(target, "crs")      <- crs
    return(as.SpatVector(target)[, 1])
  }

  .projectData <- function(x, target) {
    if (inherits(x, "SpatVector")) {
      elementSize <- terra::expanse(x, unit = "ha")
      for (i in which(terra::datatype(x) == "double")) {
        x[[i]] <- x[[i]] / elementSize
      }
      out <- list()
      for (i in seq_len(dim(target)[1])) {
        out[[i]] <- terra::intersect(target[i, ], x)
      }
      out <- do.call(rbind, out)
      elementSize <- terra::expanse(out, unit = "ha")
      for (i in which(terra::datatype(out) == "double")) {
        out[[i]] <- out[[i]] * elementSize
      }
      out <- terra::aggregate(out, by = ".id", fun = "sum", count = FALSE, na.rm = TRUE)
      names(out) <- sub("^sum\\_", "", names(out))
    } else if (inherits(x, "SpatRaster")) {
      out  <- terra::extract(x, target, "sum", bind = TRUE, na.rm = TRUE)
    }
    return(out)
  }

  .remap <- function(x, map) {
    # reduce categories to minimum based on supplied mappings
    map <- map[map$reference %in% getItems(x, dim = 3), ]
    x <- toolAggregate(x, map, from = "reference", to = "merge", dim = 3)
    return(x)
  }

  .getToolboxSpatRaster <- function(map) {
    x <- read.magpie(system.file("extdata/toolbox2000.mz", package = "mrdownscale"))
    # reduce categories to minimum based on supplied mappings
    x <- .remap(x, map)
    x <- as.SpatRaster(x)
    return(x)
  }

  .getLUH2SpatRaster <- function(map) {
    luh2 <- read.magpie(system.file("extdata/luh2015.mz", package = "mrdownscale"))
    # reduce categories to minimum based on supplied mappings
    luh2 <- .remap(luh2, map)
    luh2 <- as.SpatRaster(luh2)
    return(luh2)
  }

  toolbox <- .getToolboxSpatRaster(map)
  luh2 <- .getLUH2SpatRaster(map)

  # project toolbox and luh2 data on x
  target <- .getTarget(geometry, crs)
  # if a low res region (e.g. magpie cluster) hits only NA cells in toolbox, the result will be NA
  ptoolbox <- .projectData(toolbox, target)
  pluh2 <- .projectData(luh2, target)

  # convert to magclass
  mluh2 <- as.magpie(pluh2, spatial = which(terra::datatype(pluh2) != "double"))
  mtoolbox  <- as.magpie(ptoolbox, spatial = which(terra::datatype(ptoolbox) != "double"))
  mtoolbox[is.na(mtoolbox)] <- 0

  .dummy <- function(x, map, availableItems) {
    # generate empty dummy for missing categories
    missing <- setdiff(map$merge, availableItems)
    dummy <- x[, , rep(1, length(missing))]
    dummy[, , ] <- 0
    getItems(dummy, dim = 3) <- missing
    return(dummy)
  }

  out <- mbind(mluh2, mtoolbox, .dummy(mtoolbox, map, c(getItems(mluh2, dim = 3),
                                                        getItems(mtoolbox, dim = 3))))
  out <- toolFixWeight(out, map, from = "dataInput", to = "merge", dim = 3)
  attr(out, "crs") <- crs
  attr(out, "geometry") <- geometry

  # check data for consistency
  toolExpectTrue(identical(unname(getSets(out)[1]), "id"), "Dimensions are named correctly")
  toolExpectTrue(setequal(getItems(out, dim = 3), map$merge), "Land categories match merged categories")

  return(list(x = out,
              isocountries = FALSE,
              unit = "ha",
              min = 0,
              description = "Weights for dissagregation inputs to reference categories"))
}
