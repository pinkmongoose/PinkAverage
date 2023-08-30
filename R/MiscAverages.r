#' @include PinkAverage-package.r


#' @rdname pseudomedian
#' @title Pseudomedian
#' @description Obtain the pseudomedian of an object. This is the median of the set midpoints of all pairwise combinations of data within the sample, these midpoints being the Walsh averages.
#' @param x An `R` object that can support basic arithmetic.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @param diag A logical evaluating to `TRUE` or `FALSE` indicating whether pairs `(i,i)` should be included.
#' @param FUNC The function to apply to the Walsh averages.
#' @details This implementation includes pairs consisting of the same data point twice.
#'
#' @return The pseudomedian of these data.
#' @examples
#' pseudo_median((1:10)^2)
#' @export
pseudomedian <- function(x,na.rm=FALSE,diag=FALSE,FUNC=median) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) x <- x[!is.na(x)]
    else return(NA)
  }
  pairs <- outer(x,x,\(i,j) (i+j)/2)
  pairs[lower.tri(pairs,diag)] <- NA
  pairs <- c(pairs)
  pairs <- pairs[!is.na(pairs)]
  FUNC(pairs)
}

#' Create Bayesian average function
#'
#' Create a function to calculate a Bayesian average using a function factory. This is a mean assuming a prior mean `m` with representative weight `C`, increasing the numerator by `Cm` and the denominator by `C`.
#' @param m The prior mean.
#' @param C The prior sample size, which does not need to be an integer.
#' @details For repeat use use the associated factory function and apply it.
#' @return A function to calculate a Bayesian average of an object.
#' @examples
#' make_Bayesian_average(8,10)
#' @export
make_Bayesian_average <- function(C,m) {
  force(C)
  force(m)
  function(x,na.rm=FALSE) (C*m + sum(x,na.rm=na.rm))/(C+ifelse(na.rm,sum(!is.na(x)),length(x)))
}

#' @rdname make_Bayesian_average
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return A Bayesian average.
#' @examples
#' 
#' Bayesian_average(c(5,6,7),8,10)
#' Bayesian_average(c(NA),1,1,na.rm=TRUE)
#' @export
Bayesian_average <- function(x,C,m,na.rm=FALSE) make_Bayesian_average(C,m)(x,na.rm=na.rm)

#' Mode
#'
#' Obtain the mode (most common value) of an object.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The mode statistic. This function will return `NA` if the data contain `NA` or `NaN` and `na.rm=FALSE`, otherwise `NA` and `NaN` values are removed prior to further computation.
#' @examples 
#' mode(c(1,1,1,2,2,3))
#' @export
mode_average <- function(x,na.rm=FALSE) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) x <- x[!is.na(x)]
    else return(NA)
  }
  xunique <- unique(x)
  xunique[which.max(tabulate(match(x, xunique)))]
}
