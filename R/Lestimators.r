#' @include PinkAverage-package.r

#' @rdname n_point_L_estimator
#' @title n-point L estimators
#' @description Create an n-point L estimator function using a function factory.
#' @param q An R object representing a list of quantiles.
#' @param u An R object representing a list of weightings of the same length as `q`.
#' @return A function to calculate an n-point L estimator.
#' @examples 
#' my_trimean <- n_point_L_estimator(c(0.25,0.5,0.75), c(0.25,0.5,0.25))
#' my_trimean((1:100)^2)
#' @family L estimators
#' @export
n_point_L_estimator <- function(q,u) {
  if ((length(q)!=length(u))) stop("Number of points must match number of weights.")
  function(x,na.rm=FALSE) sum(quantile(x,q,names=FALSE,na.rm=na.rm)*u)
}

#' Create midsummary function
#'
#' Create a midsummary function using a function factory.
#' @param p The fraction (0 to 0.5) of observations to be trimmed from each end of the data before the mean is computed.
#' @return A function to calculate the midsummary on an object, produced by `n_point_L_estimator`.
#' @examples 
#' my_midsummary <- make_midsummary(0.1)
#' my_midsummary((1:10)^2)
#' @details For repeat use use associated factory function and apply it, rather than the single-function-call equivalent.
#' @family L estimators
#' @export
make_midsummary <- function(p) n_point_L_estimator(c(p,1-p),c(0.5,0.5))

#' @rdname make_midsummary
#' @param x An `R` object that can support the `quantile` function.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The midsummary statistic.
#' @examples
#' 
#' midsummary((1:10)^2,0.1)
#' @export
midsummary <- function(x,p,na.rm=FALSE) make_midsummary(p)(x,na.rm=na.rm)

#' @rdname make_midsummary
#' @return The midhinge statistic.
#' @examples 
#' midhinge((1:10)^2)
#' @export
midhinge <- make_midsummary(0.25)

#' @rdname make_midsummary
#' @return The midrange statistic.
#' @examples 
#' midrange((1:10)^2)
#' @export
midrange <- make_midsummary(0)

#' Trimean
#'
#' Obtain the trimean of an object.
#' @param x An `R` object that can support the `quantile` function.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The trimean statistic.
#' @examples 
#' trimean((1:10)^2)
#' @family L estimators
#' @export
trimean <- n_point_L_estimator(c(0.25,0.5,0.75),c(0.25,0.5,0.25))

#' @rdname trimean
#' @return The Swanson's rule statistic.
#' @examples 
#' Swansons_rule((1:10)^2)
#' @export
Swansons_rule <- n_point_L_estimator(c(0.1,0.5,0.9),c(0.3,0.4,0.3))

#' @rdname make_trimmed_range
#' @title Trimmed range
#' @description Calculates a trimmed range function using a function factory.
#' @param p The fraction (0 to 0.5) of observations to be trimmed from each end of the data before the range is computed.
#' @return A function to calculate the trimmed range on an object
#' @examples 
#' iqr <- make_trimmed_range(0.25)
#' iqr(1:100)
#' @details For repeat use use associated factory function and apply it, rather than the single-function-call equivalent.
#' @family L estimators
#' @family truncated means
#' @export
make_trimmed_range <- function(p) {
  function(x,na.rm=FALSE) {
    cna <- sum(is.na(x))
    if (cna>0) {
      if (na.rm) x <- x[!is.na(x)]
      else return(NA)
    }
    quantile(x,1-p,names=FALSE)-quantile(x,p,names=FALSE)
  }
}

#' @rdname make_trimmed_range
#' @param x An `R` object that can support the `quantile` function.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The trimmed range.
#' @examples
#' 
#' trimmed_range(1:100,0.1)
#' @export
trimmed_range <- function(x,p,na.rm=FALSE) make_trimmed_range(p)(x,na.rm=na.rm)

#' @rdname make_trimmed_range
#' @return The interdecile range.
#' @examples 
#' 
#' interdecile_range(1:100)
#' @export
interdecile_range <- make_trimmed_range(0.1)

#' @rdname make_trimmed_range
#' @return The interquartile range.
#' @examples 
#' 
#' interquartile_range(1:100)
#' @export
interquartile_range <- make_trimmed_range(0.25)

#' @rdname make_average_absolute_deviation
#' @title Average absolute deviation
#' @description Calculates an average absolute deviation function using a function factory.
#' @param CENT A function to determine a centre point of the data.
#' @param FUNC A function to apply to the deviations from the centre point of the data.
#' @param const A normalising constant to provide an estimate of the sample standard deviation.
#' @return A function to calculate the average absolute deivation on an object
#' @examples 
#' MAD <- make_average_absolute_deviation(mean,mean,sqrt(pi/2))
#' MAD(1:100)
#' @export
make_average_absolute_deviation <- function(CENT,FUNC,const=1) {
  function(x,na.rm=FALSE) {
    cna <- sum(is.na(x))
    if (cna>0) {
      if (na.rm) x <- x[!is.na(x)]
      else return(NA)
    }
    cent <- CENT(x)
    FUNC(abs(x-cent))*const
  }
}

#' @rdname make_average_absolute_deviation
#' @param x An `R` object that can support the `quantile` function.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The mean absolute deviation (mean of absolute deviations from mean, scaled by `sqrt(pi/2)`)
#' @examples
#' 
#' mean_absolute_deviation(1:100)
#' @export
mean_absolute_deviation <- function(x,na.rm=FALSE) make_average_absolute_deviation(mean,mean,sqrt(pi/2))(x,na.rm)

#' @rdname make_average_absolute_deviation
#' @return The median absolute deviation (mean of absolute deviations from median, scaled by `1.4826`)
#' @examples
#' 
#' median_absolute_deviation(1:100)
#' @export
median_absolute_deviation <- function(x,na.rm=FALSE) make_average_absolute_deviation(median,median,1.4826)(x,na.rm)

