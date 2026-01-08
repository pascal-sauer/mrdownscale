#' calcResolutionMapping
#'
#' Calculate a complete mapping from low (input dataset, clusters/countries/regions)
#' to high resolution (target dataset, grid). As a basis the mapping from the low
#' resolution clusters/countries/regions to grid cells is used.
#' Cells which are present in that mapping, but not in the target
#' dataset are discarded. Cells which are present in the target dataset, but not in the mapping are
#' added using a nearest-neighbor approach: These cells are mapped to the same low resolution
#' cluster/country/region as the closest cell which is already present in the mapping.
#'
#' @inheritParams calcLandInput
#' @param target character, the target dataset
#' @return a list including a data.frame with columns x, y, lowRes, country
#'
#' @author Pascal Sauer
calcResolutionMapping <- function(input, target) {
  if (target == "luh3pltns") {
    targetGrid <- readSource("LUH3", subtype = "states", subset = 2000)
  } else {
    targetGrid <- calcOutput("LandTargetComplete", target = target, aggregate = FALSE)
  }

  if (input == "magpie") {
    clustermap <- readSource("MagpieFulldataGdx", subtype = "clustermap")
    coords <- strsplit(clustermap$cell, "\\.")
    xCoords <- vapply(coords, function(x) as.double(sub("p", ".", x[1])), double(1))
    yCoords <- vapply(coords, function(x) as.double(sub("p", ".", x[2])), double(1))
    clustermap$cellOriginal <- sub("\\.[A-Z]{3}$", "", clustermap$cell)
    clustermap <- cbind(x = xCoords, y = yCoords, clustermap[, -which(colnames(clustermap) == "cell")])
    colnames(clustermap)[colnames(clustermap) == "cluster"] <- "lowRes"

    mapping <- toolResolutionMapping(clustermap, targetGrid)

    toolExpectTrue(all(mapping$cellOriginal %in% clustermap$cellOriginal),
                   "a subset of input cells is mapped")
  } else if (input == "witch") {
    mapping <- readSource("WITCH", subtype = "resolutionMapping")
    mapping <- mapping[, setdiff(colnames(mapping), "witch17")]
  } else if (input == "coffee") {
    mapping <- readSource("COFFEE", subtype = "regionMapping", convert = FALSE)

    magpie <- calcOutput("ResolutionMapping", input = "magpie", target = target, aggregate = FALSE)
    magpie <- magpie[, setdiff(colnames(magpie), c("region", "lowRes"))]

    mapping <- merge(magpie, mapping, by = "country")
  } else {
    stop("Unsupported input type \"", input, "\"")
  }

  toolExpectTrue(setequal(colnames(mapping), c("x", "y", "lowRes", "region", "country",
                                               "global", "cellOriginal", "cell")),
                 "resolution mapping has the expected columns")
  coords <- terra::crds(targetGrid, df = TRUE)
  allTargetCells <- paste0(sub("\\.", "p", coords$x),
                           ".",
                           sub("\\.", "p", coords$y))
  toolExpectTrue(setequal(mapping$cell, allTargetCells),
                 "all target cells are mapped")

  return(list(x = mapping,
              class = "data.frame",
              unit = NA,
              description = "mapping of high to low resolution and country"))
}
