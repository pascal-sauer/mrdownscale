# TODO documentation
calcNonlandLowRes <- function(input = "magpie", target = "luh2mod") {
  xInput <- calcOutput("NonlandHarmonizedCategories", input = input, aggregate = FALSE)

  # get target data in spatial resolution of input data
  xTarget <- calcOutput("NonlandTarget", target = target, aggregate = FALSE)
  ref <- as.SpatVector(xInput[, 1, 1])[, c(".region", ".id")]
  xTarget <- terra::extract(xTarget, ref, sum, na.rm = TRUE, bind = TRUE)
  xTarget <- as.magpie(xTarget)

  stopifnot(setequal(getItems(xInput, 3), getItems(xTarget, 3)))
  xTarget <- xTarget[, , getItems(xInput, 3)] # harmonize order of dim 3

  # TODO checks

  return(list(x = xTarget,
              isocountries = FALSE,
              unit = "Mha",
              min = 0,
              description = "Land target data at the same low resolution as the input dataset for harmonization"))
}
