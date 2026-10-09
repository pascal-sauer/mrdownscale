#' plotHarmonizedNonLand
#'
#' Plot harmonized nonland data comparisons against each other and the raw input,
#' per group, as global totals per category and harmonization. Harmonized
#' fertilizer rates in kg ha-1 yr-1 are converted back to Tg yr-1 using the
#' matching uncleaned harmonized land area. PNGs are written to the current
#' working directory.
#'
#' @param harmonizations character vector of harmonizations to plot, any of
#' "absoluteChanges", "inputUnharmonized", "fadeForest", "fade"
#' @param input name of the nonland input source
#' @param target name of the nonland target source
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets.
#' For harmonization = "absoluteChanges" both values must be set to the same
#' harmonization year, which must match a time step in both datasets.
#' @return Invisibly, the paths of the written PNG files (in the working dir).
#' @author Pascal Sauer
#' @export
plotHarmonizedNonLand <- function(harmonizations = c("absoluteChanges", "inputUnharmonized"),
                                  input = "magpie", target = "luh3",
                                  harmonizationPeriod = c(2020, 2050)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" is needed for plotHarmonizedNonLand(). Please install it.",
         call. = FALSE)
  }
  if (length(harmonizationPeriod) != 2 || anyNA(harmonizationPeriod) ||
        any(round(harmonizationPeriod) != harmonizationPeriod)) {
    stop("harmonizationPeriod must be two integer values")
  }
  stopifnot(all(harmonizations %in% c("absoluteChanges", "inputUnharmonized", "fadeForest", "fade")))

  clean <- function(x) {
    x[abs(x) < 1e-12] <- 0
    x
  }

  xInput <- clean(calcOutput("NonlandInputRecategorized", input = input, target = target, aggregate = FALSE))

  # full harmonization chains, cached
  # the harmonized land data is needed to convert the harmonized fertilizer
  # rates in kg ha-1 yr-1 back to Tg yr-1
  methods <- list()
  fertilizerTg <- list()
  for (harmonization in harmonizations) {
    if (harmonization == "inputUnharmonized") {
      methods[[harmonization]] <- xInput
      fertilizerTg[[harmonization]] <- xInput[, , "fertilizer"]
    } else {
      # calcLandHarmonized requires two equal values for absoluteChanges
      period <- if (harmonization == "absoluteChanges") harmonizationPeriod[c(1, 1)] else harmonizationPeriod
      methods[[harmonization]] <- clean(suppressMessages(calcOutput("NonlandHarmonized", input = input,
                                                                    target = target,
                                                                    harmonizationPeriod = period,
                                                                    harmonization = harmonization,
                                                                    aggregate = FALSE)))
      # deliberately uncleaned, it is used as denominator for the fertilizer conversion
      landHarmonized <- suppressMessages(calcOutput("LandHarmonized", input = input, target = target,
                                                    harmonizationPeriod = period,
                                                    harmonization = harmonization,
                                                    aggregate = FALSE))
      fertilizerTg[[harmonization]] <- clean(toolFertilizerTg(methods[[harmonization]][, , "fertilizer"],
                                                              landHarmonized))
    }
  }

  items <- getItems(xInput, 3)
  forest <- c("primf", "secmf", "secyf", "pltns")
  other <- c("primn", "secnf")
  groups <- list(
    woodHarvestArea = paste0("wood_harvest_area.", forest),
    woodHarvestOther = paste0("wood_harvest_area.", other),
    bioh = paste0("bioh.", forest),
    biohOther = paste0("bioh.", other),
    harvestWeightType = grep("^harvest_weight_type\\.", items, value = TRUE),
    fertilizer = grep("^fertilizer\\.", items, value = TRUE)
  )
  stopifnot(setequal(unlist(groups, use.names = FALSE), items))

  areaInfo <- list(scale = 1, ylab = "wood harvest area [Mha yr-1]")
  biohInfo <- list(scale = 10^-9, ylab = "wood harvest [Tg C yr-1]")
  groupInfo <- list(
    woodHarvestArea = areaInfo,
    woodHarvestOther = areaInfo,
    bioh = biohInfo,
    biohOther = biohInfo,
    harvestWeightType = list(scale = 10^-9, ylab = "harvest weight [Tg C yr-1]"),
    fertilizer = list(scale = 1, ylab = "fertilizer [Tg yr-1]")
  )

  # category labels are shown without their group prefix,
  # secondary young and mature forest are shown together as secdf
  stripPrefix <- function(item) sub("^[^.]+\\.", "", item)
  mergeSecdf <- function(item) sub("^sec(mf|yf)$", "secdf", stripPrefix(item))
  itemMaps <- list(woodHarvestArea = mergeSecdf, bioh = mergeSecdf)

  # sum magpie over cells and categories -> data.frame Year, Item, Value
  globalYearSums <- function(x) {
    d <- as.data.frame(dimSums(dimSums(x, dim = 1), dim = 1))
    setNames(data.frame(Year = as.integer(as.character(d$Year)),
                        Item = paste(d$Data1, d$Data2, sep = "."),
                        Value = d$Value), c("Year", "Item", "Value"))
  }

  files <- character(0)
  for (name in names(groups)) {
    map <- if (is.null(itemMaps[[name]])) stripPrefix else itemMaps[[name]]
    d <- do.call(rbind, lapply(names(methods), function(method) {
      x <- if (name == "fertilizer") fertilizerTg[[method]] else methods[[method]][, , groups[[name]]]
      d <- globalYearSums(x * groupInfo[[name]]$scale)
      d$Item <- map(d$Item)
      d$Method <- method
      d
    }))
    files <- c(files, plotHarmonizationGroup("nonland", name, d, ylab = groupInfo[[name]]$ylab,
                                             harmonizationPeriod = harmonizationPeriod))
  }
  return(invisible(files))
}
