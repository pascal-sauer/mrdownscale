#' toolDownscaleAvlCropland
#'
#' Downscaling method using \code{luscale::interpolateAvlCroplandWeighted} as
#' downscaling function.
#'
#' @param x magclass containing land to be downscaled
#' @param xTarget magclass target land use dataset for initialization year
#' @param xTargetLowRes magclass target land use dataset for initialization year in low resolution (like x)
#' @param mapping mapping between spatial resolutions of \code{x} and \code{xTarget}
#' @return downscaled land use dataset
#' @author Pascal Sauer
toolDownscaleAvlCropland <- function(x, xTarget, xTargetLowRes, mapping) {
  # TODO documentation, caveats
  mapping$cluster <- mapping$lowRes
  mapping <- mapping[, c("cell", "cluster")]

  stopifnot(setequal(mapping$cluster, getItems(x, 1)),
            setequal(mapping$cell, getItems(xTarget, 1)),
            setequal(getItems(x, 3), getItems(xTarget, 3)),
            nyears(xTarget) == 1)

  crops <- grep("^c[34](ann|per|nfx)", getItems(xTarget, 3), value = TRUE)

  # TODO get this on 0.25 deg via readZabel, calcAvlcropland
  avl_cropland_hr <- magclass::setNames(dimSums(xTarget[, , crops]) * 1.2, "q33_marginal")

  variableMap <- rbind(data.frame(fine = crops, magpie = "crop"),
                       data.frame(fine = c("primn", "secdn"), magpie = "other"),
                       data.frame(fine = c("pastr", "range"), magpie = "past"),
                       data.frame(fine = "primf", magpie = "primforest"),
                       data.frame(fine = c("secdf", "pltns_excl_added_treecover", "pltns_added_treecover"), magpie = "secdforest"))
  variableMap <- rbind(variableMap,
                       data.frame(fine = setdiff(getItems(x, 3), variableMap$fine),
                                  magpie = setdiff(getItems(x, 3), variableMap$fine)))

  remapVars <- function(a, weight = NULL) {
    return(toolAggregate(a, variableMap, weight = weight, dim = 3))
  }

  landTargetLowRes <- remapVars(landTargetLowRes[, getYears(xTarget), getItems(xTarget, 3)])

  "!# @monitor luscale::interpolateAvlCroplandWeighted"
  out <- luscale::interpolateAvlCroplandWeighted(x = remapVars(x)[, -1, ],
                                                 x_ini_lr = landTargetLowRes,
                                                 x_ini_hr = remapVars(xTarget),
                                                 avl_cropland_hr = avl_cropland_hr,
                                                 map = mapping,
                                                 urban_land_hr = "static",
                                                 marginal_land = "q33_marginal",
                                                 land_consv_hr = NULL, # TODO
                                                 snv_pol_shr = 0, # TODO check as in https://github.com/magpiemodel/magpie/blob/f5c2d9dead80a30ddcd0ed8593390dddbb35e576/scripts/output/extra/disaggregation.R#L296
                                                 year_ini = getYears(xTarget))

  summary(remapVars(xTarget) - out[, 1, ]) # diff of 0.077277...

  out <- remapVars(out, xTarget + 10^-10) # TODO clean up 10^-10

  # land before and after interpolateAvlCroplandWeighted does not match, this is a known problem and needs to be fixed
  # in interpolateAvlCroplandWeighted, scale to make it fit
  scalingFactor <- dimSums(x, 1) / dimSums(out, 1)
  if (!all(is.finite(scalingFactor))) {
    browser()
  }
  out <- out * scalingFactor
  # TODO this scaling does not solve the problem yet
  getSets(out)[1:2] <- c("x", "y")
  return(out)
}
