toolCheckFertilizer <- function(fertilizer, landMha = NULL, threshold = 1200) {
  # assume fertilizer is in kg ha-1 yr-1 if landMha is NULL, in Tg yr-1 otherwise
  if (!is.null(landMha)) {
    fertilizer <- toolFertilizerKgPerHa(fertilizer, landMha)
  }
  toolExpectTrue(max(fertilizer) <= threshold,
                 paste0("Fertilizer application is <= ", threshold, " kg ha-1 yr-1 (max: ",
                        signif(max(fertilizer), 3), ")"))
}
