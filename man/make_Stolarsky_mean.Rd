% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Bivariate.r
\name{make_Stolarsky_mean}
\alias{make_Stolarsky_mean}
\alias{Stolarsky_mean}
\title{Stolarsky mean}
\usage{
make_Stolarsky_mean(p)

Stolarsky_mean(x, y, p)
}
\arguments{
\item{p}{The power for the mean.}

\item{x, y}{The numbers to calculate the mean of.}
}
\value{
A function with two parameters \verb{(x,y)} to calculate a Stolarsky mean, except for the following cases. For \code{p==-Inf} the \code{min} function is returned. For \code{p==-1} the \code{geometric_mean2} function is returned. For \code{p==0} the \code{logarithmic_mean} function is returned. For \code{p==1} the \code{identric_mean} function is returned. For \code{p==2} the \code{arithmetic_mean2} function is returned. For \code{p==Inf} the \code{max} function is returned.

A Stolarsky mean.
}
\description{
Create a Stolarsky mean function to apply to two numbers, using a function factory. For a power \code{p}, the Stolarsky mean interpolates between certain other means such as the geometric mean (\code{p==-1}) and the arithmetic mean (\code{p==2}) for two numbers.
}
\details{
For repeat use use associated factory function and apply it, rather than the single-function-call equivalent.
}
\examples{
sm <- make_Stolarsky_mean(2.5)
sm(1,2)

Stolarsky_mean(1,2,2.5)
}
\references{
Stolarsky, Kenneth B. (1975). Generalizations of the logarithmic mean. Mathematics Magazine, 48: 87–92. doi:10.2307/2689825.
}
\seealso{
Other bivariate means: 
\code{\link{Heronian_mean}()},
\code{\link{arithmetic_geometric_mean}()},
\code{\link{arithmetic_mean2}()},
\code{\link{identric_mean}()},
\code{\link{logarithmic_mean}()},
\code{\link{make_Heinz_mean}()},
\code{\link{make_Seiffert_like_mean}()},
\code{\link{make_iterated_mean}()}
}
\concept{bivariate means}
