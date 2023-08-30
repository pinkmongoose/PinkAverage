#' @include PinkAverage-package.r

#' @rdname make_truncated_mean
#' @title Create truncated or Winsorized mean function
#' @description Create a mean function allowing for either truncation or Winsorizing of data, including interpolation, using a function factory. For a truncated mean, a size of tail is specified and data are excluded. For a Winsorized mean, a size of tail is specified and data in the tails are replaced by the nearest data point. Where interpolation is specified, means are interpolated between the nearest possible data point to the specified quantile, counting inclusively and exclusively.
#' @param FUNC The function to apply to the truncated or Winsorized data, usually `mean`.
#' @param tail The size of tail at each end of the sorted data to be truncated or Winsorized.
#' @param type The type of tail. Can be `count`, `quartile` (0-4), `quantile` (0-1), `decile` (0-10), or `percentile` (0-100).
#' @param Winsorize The tails are Winsorized if true, otherwise excluded.
#' @param interpolate Tails are conservative/exclusive if false, otherwise interpolation is performed.
#' @return A function with parameters `(x, na.rm=FALSE)` to calculate a truncated or Winsorized mean on an object. These functions will return `NA` if the data contain `NA` or `NaN` and `na.rm=FALSE`, otherwise `NA` and `NaN` values are removed prior to further computation.
#' @examples 
#' iqmean <- make_truncated_mean(mean, 1, type="quartile", interpolate=TRUE)
#' iqmean((1:10)^2)
#' @family truncated means
#' @export
make_truncated_mean <- function(FUNC, tail, type="count", Winsorize=FALSE, interpolate=TRUE) function(x,na.rm=FALSE) {
  cna <- sum(is.na(x))
  if (cna>0) {
    if (na.rm) x <- x[!is.na(x)]
    else return(NA)
  }
  if (type=="quartile") tail <- tail*0.25*length(x)
  else if (type=="quantile") tail <- tail*length(x)
  else if (type=="decile") tail<- tail*length(x)/10
  else if (type=="percentile") tail <- tail*length(x)/100
  else if (type!="count") stop(paste("Truncation type '",type,"' not recognised."))
  if (!is.vector(x)) stop("x is not a vector in truncated_mean.")
  if (interpolate) {
    if (ceiling(tail)+1>floor(length(x)-tail)) stop("Can't truncate more than whole dataset in truncated_mean.")
  } else {
    if (floor(tail)+1>ceiling(length(x)-tail)) stop("Can't truncate more than whole dataset in truncated_mean.")
  }
  x <- sort(x)
  if (Winsorize) x1 <- sapply(seq_along(x),
                              \(i) ifelse(i>=floor(tail)+1, ifelse(i<=ceiling(length(x)-tail), x[i], x[ceiling(length(x)-tail)]), x[floor(tail)+1]))
  else x1 <- x[(floor(tail)+1):(ceiling(length(x)-tail))]
  #print(x1)
  if (interpolate) {
    if (Winsorize) x2 <- sapply(seq_along(x),
                                \(i) ifelse(i>=ceiling(tail)+1, ifelse(i<=floor(length(x)-tail), x[i], x[floor(length(x)-tail)]), x[ceiling(tail)+1]))
    else x2 <- x[(ceiling(tail)+1):(floor(length(x)-tail))]
    #print(x2)
    frac <- tail%%1
    return(FUNC(x1)*(1-frac) + FUNC(x2)*frac)
  } else {
    return(FUNC(x1))
  }
}

#' @rdname modified_mean
#' @title Modified mean
#' @description Obtain the modified mean of an an object, also known as the Olympic average. This function is derived from the function factory `make_truncated_mean`. The modified mean excludes the single largest and smallest observations.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The modified mean.
#' @examples 
#' modified_mean((1:10)^2)
#' @family truncated means
#' @export
modified_mean <- make_truncated_mean(mean, 1, "count", FALSE, FALSE)

#' @rdname modified_mean
#' @export
Olympic_average <- modified_mean

#' Truncated mean
#'
#' Obtain the truncated mean of an an object, which excludes the tails of the data. This function is derived from the function factory `make_truncated_mean`. Interpolation is specified where the length of `x` is not a multiple of 4.
#' @param x An `R` object.
#' @param tail The tails of the distribution to be removed, specified as a quantile (0-1).
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The truncated mean.
#' @examples 
#' truncated_mean((1:10)^2,0.1)
#' @family truncated means
#' @export
truncated_mean <- function(x,tail,na.rm=FALSE) make_truncated_mean(mean, tail, "quantile", FALSE, TRUE)(x,na.rm=na.rm)

#' Winsorized mean
#'
#' Obtain the Winsorized mean of an an object, which replaces the tails of the data with the nearest value outside the tails. This function is derived from the function factory `make_truncated_mean`. Interpolation is specified where the length of `x` is not a multiple of 4.
#' @param x An `R` object.
#' @param tail The tails of the distribution to be Winsorized, specified as a quantile (0-1).
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The Winsorized mean.
#' @examples 
#' Winsorized_mean((1:10)^2,0.1)
#' @family truncated means
#' @export
Winsorized_mean <- function(x,tail,na.rm=FALSE) make_truncated_mean(mean, tail, "quantile", TRUE, TRUE)(x,na.rm=na.rm)

#' @rdname Winsorized_mean
#' @examples Winsorised_mean((1:10)^2,0.1)
#' @export
Winsorised_mean <- Winsorized_mean

#' Interquartile mean
#'
#' Obtain the interquartile mean of an an object, which excludes the first and fourth quartiles of data. This function is derived from the function factory `make_truncated_mean`. Interpolation is specified where the length of `x` is not a multiple of 4.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The interquartile mean.
#' @examples 
#' interquartile_mean((1:10)^2)
#' @family truncated means
#' @export
interquartile_mean <- function(x,na.rm=FALSE) make_truncated_mean(mean, 0.25, "quantile", FALSE, TRUE)(x,na.rm=na.rm)
