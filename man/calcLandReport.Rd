% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcLandReport.R
\name{calcLandReport}
\alias{calcLandReport}
\title{calcLandReport}
\usage{
calcLandReport(
  outputFormat,
  harmonizationPeriod,
  yearsSubset,
  harmonization,
  downscaling
)
}
\arguments{
\item{outputFormat}{format in which the outputs should be prepared. Options: ESM, ScenarioMIP}

\item{harmonizationPeriod}{Two integer values, before the first given
year the target dataset is used, after the second given year the input
dataset is used, in between harmonize between the two datasets}

\item{yearsSubset}{vector of years to keep in the output dataset}

\item{harmonization}{name of harmonization method, see \code{\link{toolGetHarmonizer}}}

\item{downscaling}{name of downscaling method, currently only "magpieClassic"}
}
\value{
land use data
}
\description{
Convert the downscaled land use data to the format required by the given project.
}
\author{
Pascal Sauer
}
