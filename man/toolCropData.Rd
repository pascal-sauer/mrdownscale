% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/toolCropData.R
\name{toolCropData}
\alias{toolCropData}
\title{toolCropData}
\usage{
toolCropData(landHighRes, cellArea)
}
\arguments{
\item{landHighRes}{high resolution land use data as magclass object}

\item{cellArea}{corresponding magclass object containing cell area in Mha}
}
\value{
crop data as magclass object with variables irrig_*, cpbf1_*,
cpbf2_*, * for all 5 crop types
}
\description{
For each crop type (c3ann, c3nfx, c3per, c4ann, c4per) calculate irrigation
share, 1st and 2nd generation biofuel shares and share of cell area
(corresponding to LUH3 variables c3ann, irrig_c3ann, cpbf1_c3ann,
cpbf2_c3ann and analogously for the other crop types).
}
\author{
Pascal Sauer
}
