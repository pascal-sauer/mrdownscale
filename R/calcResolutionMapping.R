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
#' @param input character, the input dataset
#' @param target character, the target dataset
#' @return a list including a data.frame with columns x, y, lowRes, countrycode
#'
#' @author Pascal Sauer
calcResolutionMapping <- function(input, target) {
  if (input == "magpie") {
    clustermap <- readSource("MagpieFulldataGdx", subtype = "clustermap")
    coords <- strsplit(clustermap$cell, "\\.")
    xCoords <- vapply(coords, function(x) as.double(sub("p", ".", x[1])), double(1))
    yCoords <- vapply(coords, function(x) as.double(sub("p", ".", x[2])), double(1))
    clustermap$cellOriginal <- sub("\\.[A-Z]{3}$", "", clustermap$cell)
    clustermap <- cbind(x = xCoords, y = yCoords, clustermap[, -which(colnames(clustermap) == "cell")])
    colnames(clustermap)[colnames(clustermap) == "cluster"] <- "lowRes"
  } else {
    stop("Unsupported input type \"", input, "\"")
  }

  targetGrid <- calcOutput("LandTarget", target = target, aggregate = FALSE)
  mapping <- toolResolutionMapping(clustermap, targetGrid)

  toolExpectTrue(setequal(colnames(mapping), c("x", "y", "lowRes", "region", "country",
                                               "global", "cellOriginal", "cell")),
                 "resolution mapping has the expected columns")
  coords <- terra::crds(targetGrid, df = TRUE)
  allTargetCells <- paste0(sub("\\.", "p", coords$x),
                           ".",
                           sub("\\.", "p", coords$y))
  toolExpectTrue(setequal(mapping$cell, allTargetCells),
                 "all target cells are mapped")
  toolExpectTrue(all(mapping$cellOriginal %in% clustermap$cellOriginal),
                 "a subset of input cells is mapped")

  return(list(x = mapping,
              class = "data.frame",
              unit = NA,
              description = "mapping of high to low resolution and countrycode"))
}
