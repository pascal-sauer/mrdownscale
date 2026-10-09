#' plotHarmonizedLand
#'
#' Plot harmonized land data comparisons against each other and the raw input,
#' per group, as global totals per category and harmonization. PNGs are written
#' to the current working directory.
#'
#' @param harmonizations character vector of harmonizations to plot, any of
#' "absoluteChanges", "inputUnharmonized", "fadeForest", "fade"
#' @param input name of the land input source
#' @param target name of the land target source
#' @param harmonizationPeriod Two integer values, before the first given
#' year the target dataset is used, after the second given year the input
#' dataset is used, in between harmonize between the two datasets.
#' For harmonization = "absoluteChanges" both values must be set to the same
#' harmonization year, which must match a time step in both datasets.
#' @return Invisibly, the paths of the written PNG files (in the working dir).
#' @author Pascal Sauer
#' @export
plotHarmonizedLand <- function(harmonizations = c("absoluteChanges", "inputUnharmonized"),
                               input = "magpie", target = "luh3",
                               harmonizationPeriod = c(2020, 2050)) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package \"ggplot2\" is needed for plotHarmonizedLand(). Please install it.",
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

  xInput <- clean(calcOutput("LandInputRecategorized", input = input, target = target, aggregate = FALSE))

  # full harmonization chains, cached
  methods <- list()
  for (harmonization in harmonizations) {
    if (harmonization == "inputUnharmonized") {
      xTarget <- calcOutput("LandTargetLowRes", input = input, target = target,
                            endOfHistory = harmonizationPeriod[1], aggregate = FALSE)
      methods[[harmonization]] <- clean(toolEqualizeArea(xInput, xTarget[, harmonizationPeriod[1], ]))
    } else {
      # calcLandHarmonized requires two equal values for absoluteChanges
      period <- if (harmonization == "absoluteChanges") harmonizationPeriod[c(1, 1)] else harmonizationPeriod
      methods[[harmonization]] <- clean(suppressMessages(calcOutput("LandHarmonized", input = input,
                                                                    target = target,
                                                                    harmonizationPeriod = period,
                                                                    harmonization = harmonization,
                                                                    aggregate = FALSE)))
    }
  }

  # urban is excluded from all plots
  items <- getItems(xInput, 3)
  groups <- list(
    forest = intersect(c("pltns", "primf", "secdf"), items),
    primnSecdn = intersect(c("primn", "secdn"), items),
    cropland = grep("rainfed|irrigated", items, value = TRUE),
    pastureRangeland = intersect(c("pastr", "range"), items)
  )
  other <- setdiff(items, c(unlist(groups, use.names = FALSE), "urban"))
  if (length(other) > 0) {
    groups$other <- other
  }

  # merges irrigated/rainfed into one entry per crop (and per biofuel generation),
  # shows all 1st gen biofuel crops together (it is zero in the magpie input)
  croplandGroup <- function(item) {
    item <- sub("_(irrigated|rainfed)(_|$)", "\\2", item)
    sub("^[a-z0-9]+_biofuel_1st_gen$", "biofuel_1st_gen", item)
  }
  itemMaps <- list(cropland = croplandGroup)

  # sum magpie over cells -> data.frame Year, Item, Value
  globalYearSums <- function(x) {
    d <- as.data.frame(dimSums(x, dim = 1))
    setNames(data.frame(Year = as.integer(as.character(d$Year)), Item = d$Data1,
                        Value = d$Value), c("Year", "Item", "Value"))
  }

  files <- character(0)
  for (name in names(groups)) {
    map <- if (is.null(itemMaps[[name]])) identity else itemMaps[[name]]
    d <- do.call(rbind, lapply(names(methods), function(method) {
      d <- globalYearSums(methods[[method]][, , groups[[name]]])
      d$Method <- method
      d
    }))
    d$Item <- map(d$Item)
    files <- c(files, plotHarmonizationGroup("land", name, d, ylab = "area [Mha]",
                                             harmonizationPeriod = harmonizationPeriod))
  }
  return(invisible(files))
}
