#' @include PinkAverage-package.r

#' @rdname circular_mean
#' @title Circular mean
#' @description Obtain the circular mean of an object.
#' @param x An `R` object that can support trigonometric functions, representing a series of angles measured in radians or number over a given modulus.
#' 
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The circular mean of these angles, between `-pi` and `+pi`.
#' @details A danger of circular mean functions is that where data points are perfectly balanced around the unit circle, the circular mean becomes singular, in which case zero or effectively random results can be returned dependent upon rounding.
#' @examples 
#' circular_mean(1:5)
#' @export
circular_mean <- function(x,na.rm=FALSE) atan2(mean(sin(x),na.rm=na.rm),mean(cos(x),na.rm=na.rm))

#' @rdname circular_mean
#' @param m The modulus.
#' @return The modular mean of these numbers, between 0 and `m`.
#' @examples 
#' 
#' modular_mean(c(2,8),10)
#' @export 
modular_mean <- function(x,m=1,na.rm=FALSE) (circular_mean(2*pi*x/m)*m/(2*pi)+m)%%m

#' @rdname circular_mean
#' @export
degrees_mean <- function(x,na.rm=FALSE) modular_mean(x,360,na.rm)

#' @rdname circular_mean
#' @export
gradians_mean <- function(x,na.rm=FALSE) modular_mean(x,400,na.rm)

#' @rdname circular_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_circular_mean(1:5,c(0,1,1,1,0))
#' @export
weighted_circular_mean <- function(x,w,na.rm=FALSE) atan2(weighted_mean(sin(xa),w,na.rm=na.rm),weighted_mean(cos(x),w,na.rm=na.rm))

#' @rdname circular_mean
#' @examples 
#' 
#' weighted_modular_mean(c(2,8,9),c(1,1,0),10)
#' @export
weighted_modular_mean <- function(x,w,m,na.rm=FALSE) (weighted_circular_mean(2*pi*x/m,w)*m/(2*pi)+m)%%m

#' @rdname circular_mean
#' @export
weighted_degrees_mean <- function(x,w,na.rm=FALSE) modular_mean(x,w,360,na.rm)

#' @rdname circular_mean
#' @export
weighted_gradians_mean <- function(x,w,na.rm=FALSE) modular_mean(x,w,400,na.rm)
