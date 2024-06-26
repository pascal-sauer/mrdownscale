#' downloadLUH2v2h
#'
#' Download the LUH2v2h dataset (states.nc, transitions.nc,
#' management.nc, staticData_quarterdeg.nc).
#'
#' @param subtype one of states, management, transitions, cellArea
#' @return metadata list with URL, DOI, title, description, author, unit, version, release date
#' @author Pascal Sauer
downloadLUH2v2h <- function(subtype = NULL) {
  fileNames <- list(states = c("states.nc", "staticData_quarterdeg.nc"),
                    transitions = c("transitions.nc", "staticData_quarterdeg.nc"),
                    management = "management.nc",
                    cellArea = "staticData_quarterdeg.nc")

  fileNames <- toolSubtypeSelect(subtype, fileNames)

  # increase timeout drastically because the whole file needs to be downloaded in that time (16GB)
  withr::local_options(timeout = max(3e6, getOption("timeout")))
  for (fileName in fileNames) {
    utils::download.file(paste0("https://luh.umd.edu/LUH2/LUH2_v2h/", fileName), fileName, mode = "wb")
  }

  return(list(url          = "https://luh.umd.edu/index.shtml",
              doi          = "10.5194",
              title        = paste("LUH2 v2h", subtype),
              description  = paste("The goal of the Land-Use Harmonization (LUH2) project is to",
                                   "prepare a harmonized set of land-use scenarios that smoothly",
                                   "connects the historical reconstructions of land-use with the",
                                   "future projections in the format required for ESMs."),
              author       = "Hurtt et al.",
              unit         = "1",
              version      = "v2h",
              release_date = "2016-10-14",
              license      = NA,
              reference    = NA))
}
