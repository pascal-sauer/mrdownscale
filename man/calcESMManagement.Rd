% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcESMManagement.R
\name{calcESMManagement}
\alias{calcESMManagement}
\title{calcESMManagement}
\usage{
calcESMManagement(harmonizationPeriod, yearsSubset)
}
\arguments{
\item{harmonizationPeriod}{Two integer values, before the first given
year the target dataset is used, after the second given year the input
dataset is used, in between harmonize between the two datasets}

\item{yearsSubset}{remove years from the returned data which are not in yearsSubset}
}
\value{
ESM compatible management information
}
\description{
Create ESM compatible management information
}
\author{
Pascal Sauer, Jan Philipp Dietrich
}
