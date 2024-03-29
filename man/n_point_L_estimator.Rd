% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Lestimators.r
\name{n_point_L_estimator}
\alias{n_point_L_estimator}
\title{n-point L estimators}
\usage{
n_point_L_estimator(q, u)
}
\arguments{
\item{q}{An R object representing a list of quantiles.}

\item{u}{An R object representing a list of weightings of the same length as \code{q}.}
}
\value{
A function to calculate an n-point L estimator.
}
\description{
Create an n-point L estimator function using a function factory.
}
\examples{
my_trimean <- n_point_L_estimator(c(0.25,0.5,0.75), c(0.25,0.5,0.25))
my_trimean((1:100)^2)
}
\seealso{
Other L estimators: 
\code{\link{make_midsummary}()},
\code{\link{make_trimmed_range}()},
\code{\link{trimean}()}
}
\concept{L estimators}
