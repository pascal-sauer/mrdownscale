% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolResolutionMapping.R
\name{toolResolutionMapping}
\alias{toolResolutionMapping}
\title{toolResolutionMapping}
\usage{
toolResolutionMapping(mapping, targetGrid)
}
\arguments{
\item{mapping}{a data.frame with columns x, y, lowRes}

\item{targetGrid}{a terra SpatRaster with the target resolution}
}
\value{
a data.frame with columns x, y, lowRes, countrycode
}
\description{
See description of \code{\link{calcResolutionMapping}}. Here we are
assuming target resolution is finer than what mapping already provides.
}
\author{
Pascal Sauer
}
