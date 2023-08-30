#' @include PinkAverage-package.r

#' @rdname make_Lehmer_mean
#' @title Lehmer mean functions
#' @description Create a Lehmer mean function using a function factory. For a power `p`, the Lehmer mean is the ratio of the sum of the data points raised to the power of `p` and the sum of the data points raised to the power of `p-1`. The mean interpolates between the arithmetic mean for `p=1` and the harmonic mean for `p=0` and includes other means at other values of `p`.
#' @param p The power used in the numerator of the expression.
#' @details For very high or low values of `p` overflow or underflow during internal calculations may generate unexpected results. For repeat use use associated factory function and apply it, rather than the single-function-call equivalent. Removal of `NA` and `NaN` values is performed by `mean`, `min`, `max` or `sum`, therefore will remove any such values introduced by manipulating the data.
#' @return A function with parameters `(x, na.rm=FALSE)` to calculate a Lehmer mean of an object, except for the following cases. For `p==-Inf` the `min` function is returned. For `p==0` the `harmonic_mean` function is returned. For `p==1` the `mean` function is returned. For `p==2` a function optimised to avoid power functions is returned. For `p==Inf` the `max` function is returned.
#' @examples 
#' contramean <- make_Lehmer_mean(2)
#' contramean(c(1,2,3))
#' @family Lehmer means
#' @references 
#' Bullen, P. S. (1987). Handbook of means and their inequalities. Springer.
#' @export
make_Lehmer_mean <- function(p) {
  force(p)
  if (p==-Inf) return(min)
  else if (p==0) return(harmonic_mean)
  else if (p==1) return(mean)
  else if (p==2) return(function(x,na.rm=FALSE) sum(x*x,na.rm=na.rm)/sum(x,na.rm=na.rm))
  else if (p==Inf) return(max)
  function(x,na.rm=FALSE) sum(x^p,na.rm=na.rm) / sum(x^(p-1),na.rm=na.rm)
}

#' @rdname make_Lehmer_mean
#' @return A function with parameters `(x, w, na.rm=FALSE)` to calculate a weighted Lehmer mean of an object, with associated weights, except for the following cases. For `p==-Inf` the `weighted_min` function is returned. For `p==0` the `weighted_harmonic_mean` function is returned. For `p==1` the `weighted_mean` function is returned. For `p==2` a function optimised to avoid power functions is returned. For `p==Inf` the `weighted_max` function is returned.
#' @examples
#' 
#' contrameanw <- make_weighted_Lehmer_mean(2)
#' contrameanw(c(1,2,3,4),c(1,1,1,0))
#' @export
make_weighted_Lehmer_mean <- function(p) {
  force(p)
  if (p==-Inf) return(weighted_min)
  else if (p==0) return(weighted_harmonic_mean)
  else if (p==1) return(weighted_mean)
  else if (p==2) return(function(x,w,na.rm=FALSE) sum(w*x*x,na.rm=na.rm)/sum(w*x,na.rm=na.rm))
  else if (p==Inf) return(weighted_max)
  function(x,na.rm=FALSE) sum(w*(x^p),na.rm=na.rm) / sum(w*(x^(p-1)),na.rm=na.rm)
}

#' @rdname make_Lehmer_mean
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The Lehmer mean.
#' @examples
#' 
#' Lehmer_mean(c(1,2,3),2)
#' @export
Lehmer_mean <- function(x,p,na.rm=FALSE) make_Lehmer_mean(p)(x,na.rm=na.rm)

#' @rdname make_Lehmer_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_Lehmer_mean(c(1,2,3,4),c(1,1,1,0),2)
#' @export
weighted_Lehmer_mean <- function(x,w,p,na.rm=FALSE) make_weighted_Lehmer_mean(p)(x,w,na.rm=na.rm)

#' Contraharmonic mean
#'
#' Obtain the contraharmonic mean of an an object, also known as the self-weighting mean. I.e. the sum of the squared data divided by the sum of the data. This function is derived from the function factory `make_Lehmer_mean` with `p=2`.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The contraharmonic mean.
#' @examples 
#' contraharmonic_mean(1:3)
#' @references 
#' Zelen, M. (1972) Length-biased sampling and biomedical problems. In Biometric Society Meeting, Dallas, Texas.
#' @family Lehmer means
#' @export
contraharmonic_mean <- make_Lehmer_mean(2)

#' @rdname contraharmonic_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_contraharmonic_mean(1:4,c(1,1,1,0))
#' @export
weighted_contraharmonic_mean <- make_weighted_Lehmer_mean(2)
