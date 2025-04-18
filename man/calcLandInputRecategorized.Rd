% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcLandInputRecategorized.R
\name{calcLandInputRecategorized}
\alias{calcLandInputRecategorized}
\title{calcLandInputRecategorized}
\usage{
calcLandInputRecategorized(input, target)
}
\arguments{
\item{input}{name of the land input source to be used}

\item{target}{name of the land target source to be used}
}
\description{
Computes the land input data in target land categories. Splitting of land
categories is performed under use of internal land weights reflecting the
prevalence of a certain land category in the given area.
}
\details{
Mapping from input to target categories is achieved via a merge of a land input
mapping to reference categories and a mapping between land target categories and
the same reference categories. Thereby a new source or new target can be supported
by supplying a map of that new input and/or target to the reference categories.
}
\author{
Jan Philipp Dietrich, Pascal Sauer
}
