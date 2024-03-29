% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Truncated.r
\name{modified_mean}
\alias{modified_mean}
\alias{Olympic_average}
\title{Modified mean}
\usage{
modified_mean(x, na.rm = FALSE)

Olympic_average(x, na.rm = FALSE)
}
\arguments{
\item{x}{An \code{R} object.}

\item{na.rm}{A logical evaluating to \code{TRUE} or \code{FALSE} indicating whether \code{NA} values should be stripped before the computation proceeds.}
}
\value{
The modified mean.
}
\description{
Obtain the modified mean of an an object, also known as the Olympic average. This function is derived from the function factory \code{make_truncated_mean}. The modified mean excludes the single largest and smallest observations.
}
\examples{
modified_mean((1:10)^2)
}
\seealso{
Other truncated means: 
\code{\link{Winsorized_mean}()},
\code{\link{interquartile_mean}()},
\code{\link{make_trimmed_range}()},
\code{\link{make_truncated_mean}()},
\code{\link{truncated_mean}()}
}
\concept{truncated means}
