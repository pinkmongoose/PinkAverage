% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bivariate.r
\name{arithmetic_mean2}
\alias{arithmetic_mean2}
\alias{geometric_mean2}
\alias{harmonic_mean2}
\alias{contraharmonic_mean2}
\title{Bivariate mean helper functions}
\usage{
arithmetic_mean2(x, y)

geometric_mean2(x, y)

harmonic_mean2(x, y)

contraharmonic_mean2(x, y)
}
\arguments{
\item{x, y}{The numbers to calculate the mean of.}
}
\value{
The arithmetic mean.

The geometric mean.

The harmonic mean.

The contraharmonic mean.
}
\description{
Helper functions for use in iterated means, to calculate simple means with two numeric parameters.
}
\examples{
arithmetic_mean2(1,2)
}
\seealso{
Other iterated means: 
\code{\link{arithmetic_geometric_mean}()},
\code{\link{make_iterated_mean}()}

Other bivariate means: 
\code{\link{Heronian_mean}()},
\code{\link{arithmetic_geometric_mean}()},
\code{\link{identric_mean}()},
\code{\link{logarithmic_mean}()},
\code{\link{make_Heinz_mean}()},
\code{\link{make_Seiffert_like_mean}()},
\code{\link{make_Stolarsky_mean}()},
\code{\link{make_iterated_mean}()}
}
\concept{bivariate means}
\concept{iterated means}
