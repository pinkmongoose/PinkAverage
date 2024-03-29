% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MiscAverages.r
\name{mode_average}
\alias{mode_average}
\title{Mode}
\usage{
mode_average(x, na.rm = FALSE)
}
\arguments{
\item{x}{An \code{R} object.}

\item{na.rm}{A logical evaluating to \code{TRUE} or \code{FALSE} indicating whether \code{NA} values should be stripped before the computation proceeds.}
}
\value{
The mode statistic. This function will return \code{NA} if the data contain \code{NA} or \code{NaN} and \code{na.rm=FALSE}, otherwise \code{NA} and \code{NaN} values are removed prior to further computation.
}
\description{
Obtain the mode (most common value) of an object.
}
\examples{
mode(c(1,1,1,2,2,3))
}
