#' @include PinkAverage-package.r

#' @rdname weighted_mean
#' @title Weighted mean function
#' @description The core weighted.mean function treats zero-weighted infinities as zero. For internal consistency within the library, `weighted_mean` instead allows NaNs to propagate, treating the result of `0*Inf` to be singular.
#' @param x An `R` object.
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The weighted mean.
#' @examples 
#' weighted_mean(c(-Inf,1,Inf),c(0,1,0))
#' weighted.mean(c(-Inf,1,Inf),c(0,1,0))
#' @family Weighted summary statistics
#' @export
weighted_mean <- function(x,w,na.rm=FALSE) {
  if (length(x)!=length(w)) stop("Data and weights have different lengths.")
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) x <- x[!is.na(x)]
    else return(NA)
  }
  sum(x*w)/sum(w)
}

#' @rdname weighted_var
#' @title weighted variance
#' @description The weighted_var function treats zero-weighted allows NaNs to propagate, treating the result of `0*Inf` to be singular.
#' @param x An `R` object.
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @param norm A logical evaluating to `TRUE` or `FALSE` indicating whether weights should be normalised to `n-1` before computation proceeds.
#' @return The weighted variance.
#' @details There are two forms of weighted variance: probability/reliability weights, and occurrences/repeat weights. By default, `weighted_variance` assumes the former and rescales the weights back to a total of `n-1`.
#' @examples 
#' weighted_var(c(1,2,3),c(0.5,0.5,0))
#' weighted_var(c(1,2,3),c(1,1,0))
#' @family Weighted summary statistics
#' @export
weighted_var <- function(x,w,na.rm=FALSE,norm=TRUE) {
  if (length(x)!=length(w)) stop("Data and weights have different lengths.")
  if (norm) w <- (length(x)-1) * w/sum(w)
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) {
      w <- w[!is.na(w)]
      x <- x[!is.na(x)]
    }  
    else return(NA)
  }
  mu <- weighted_mean(x,w)
  sum(w*(x-mu)*(x-mu))/(sum(w)-1)
}

#' @rdname weighted_var
#' @return The weighted standard deviation.
#' @examples 
#' 
#' weighted_sd(c(1,2,3),c(0.5,0.5,0))
#' @export
weighted_sd <- function(x,w,na.rm=FALSE,norm=TRUE) sqrt(weighted_var(x,w,na.rm,norm))

#' @rdname weighted_min
#' @title Weighted minimum and maximum
#' @description
#' Generates the minimum or maximum of an object, of those values where the weightings provided by a second object are non zero. These functions are returned as an asymptotic case for weighted mean factory functions where `min` and `max` would be returned by the non-weighted factory functions.
#' @param x An `R` object.
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The minimum value, or `NaN` where this has zero weighting.
#' @details For a mean such as the Holder mean with parameter \eqn{p \to  -\infty}, the mean approaches the minimum (and vice versa for the maximum) value of the dataset for all non-zero weightings. For other means, the lowest value of the dataset does not contribute to the mean where its waiting zero. The current function treats this situation as a singularity and returns `NaN` where the minimum value has zero weights (and vice versa for the maximum).
#' @examples 
#' weighted_min(1:4,c(0,1,1,1))
#' @family Kolmogorov means
#' @export
weighted_min <- function(x,w,na.rm=FALSE) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) {
      w <- w[!is.na(x)]
      x <- x[!is.na(x)]
    } else return(NA)
  }
  xmin <- which.min(x)
  ifelse(w[xmin]!=0, x[xmin], NaN)
}

#' @rdname weighted_min
#' @return The maximum value, or `NaN` where this has zero weighting.
#' @examples 
#' 
#' weighted_max(1:4,c(1,1,1,0))
#' @export
weighted_max <- function(x,w,na.rm=FALSE) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) {
      w <- w[!is.na(x)]
      x <- x[!is.na(x)]
    } else return(NA)
  }
  xmax <- which.max(x)
  ifelse(w[xmax]!=0, x[xmax], NaN)
}
