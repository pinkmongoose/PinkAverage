#' @include PinkAverage-package.r

#' @rdname Frechet_mean
#' @title Fréchet mean functions
#' @description Also known as a Karcher mean. Calculate a mean by minimising the sum of a distance metric from the mean to the data. For different choices of distance metric this generalizes means such as the arithmetic mean and the median.
#' @param DIST a function to scale the measure \eqn{x - \bar x}, which is required to be vectorised.
#' @param ... Other parameters to be passed to `nlm`.
#' @return A function with parameters `(x, na.rm=FALSE)` to calculate a Fréchet mean of an object. These functions will return `NA` if the data contain `NA` or `NaN` and `na.rm=FALSE`, otherwise `NA` and `NaN` values are removed prior to further computation.
#' @details For repeat use use associated factory function and apply it, rather than the single-function-call equivalent.
#' @examples
#' x <- c(1,3,6,10,15) 
#' my_average_mean <- make_Frechet_mean(\(z) z*z)
#' my_average_mean(x)
#' 
#' my_average_median <- make_Frechet_mean(abs)
#' my_average_median(x)
#' @family Frechet means
#' @export
make_Frechet_mean <- function(DIST, ...) function(x,na.rm=FALSE) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) x <- x[!is.na(x)]
    else return(NA)
  }
  f <- function(p, x) sum(DIST(p-x))
  nlm(f, mean(range(x)), x, ...)$estimate
}

#' @rdname Frechet_mean
#' @return A function with parameters `(x, w, na.rm=FALSE)` to calculate a weighted Frechet mean of an object.
#' @export
make_weighted_Frechet_mean <- function(DIST, ...) function(x,w,na.rm=FALSE) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) {
      w <- w[!is.na(x)]
      x <- x[!is.na(x)]
    } else return(NA)
  }
  f <- function(p, x, w) sum(w*DIST(p-x))
  nlm(f, mean(range(x)), x, w, ...)$estimate
}

#' @rdname Frechet_mean
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return A Frechet mean.
#' @examples 
#' 
#' Frechet_mean(c(1,3,6,10,15),abs)
Frechet_mean <- function(x,DIST,na.rm=FALSE,...) make_Frechet_mean(DIST,...)(x,na.rm=na.rm)

#' @rdname Frechet_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples 
#' 
#' weighted_Frechet_mean(c(1,3,6,10,15),c(1,1,1,0,0),abs)
weighted_Frechet_mean <- function(x,w,DIST,na.rm=FALSE,...) make_weighted_Frechet_mean(DIST,...)(x,w,na.rm=na.rm)
