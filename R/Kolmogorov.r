#' @include PinkAverage-package.r
#' @include SummaryStats.r

#' @rdname make_Kolmogorov_mean
#' @title Kolmogorov mean functions
#' @description Create a quasi-arithmetic mean function using a function factory. Also known as a Kolmogorov mean or generalised _f_-mean. A summary statistic is calculated using the generic `mean` or other supplied function `FUNC` on transformed data using the function `TRA` and back transformed using the function `INV`, which are required to be vectorised.
#' @param TRA The transformation function.
#' @param INV The inverse function corresponding to TRA.
#' @param FUNC The actual function to apply to the transformed data, able to support `na.rm`.
#' @return A function with parameters `(x, na.rm=FALSE)` to calculate a Kolmogorov mean of an object.
#' @details For repeat use use associated factory function and apply it, rather than the single-function-call equivalent. Removal of `NA` and `NaN` values is performed by `mean`, therefore will remove any such values introduced by `TRA`.
#' @examples 
#' geomean <- make_Kolmogorov_mean(log,exp)
#' geomean(c(1,10,100))
#' geomean(c(1,10,100,NA), na.rm=TRUE)
#' @family Kolmogorov means
#' @export
make_Kolmogorov_mean <- function(TRA,INV,FUNC=mean) {
  force(TRA)
  force(INV)
  force(FUNC)
  function(x,na.rm=FALSE) INV(FUNC(TRA(x),na.rm=na.rm))
}

#' @rdname make_Kolmogorov_mean
#' @return A function with parameters `(x, w, na.rm=FALSE)` to calculate a weighted Kolmogorov mean of an object.
#' @examples
#'
#' geomeanw <- make_weighted_Kolmogorov_mean(log,exp)
#' geomeanw(c(1,10,100),c(0,1,1))
#' @export
make_weighted_Kolmogorov_mean <- function(TRA,INV,FUNC=weighted_mean) {
  force(TRA)
  force(INV)
  force(FUNC)
  function(x,w,na.rm=FALSE) INV(FUNC(TRA(x),w,na.rm=na.rm))
}

#' @rdname make_Kolmogorov_mean
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return A Kolmogorov mean of an object.
#' @examples
#' 
#' Kolmogorov_mean(c(1,10,100),log,exp)
#' @export
Kolmogorov_mean <- function(x,TRA,INV,na.rm=FALSE) make_Kolmogorov_mean(TRA,INV)(x,na.rm=na.rm)


#' @rdname make_Kolmogorov_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples 
#' 
#' weighted_Kolmogorov mean(c(1,10,100),c(0,1,1),log,exp)
#' @export
weighted_Kolmogorov_mean <- function(x,w,TRA,INV,na.rm=FALSE) make_weighted_Kolmogorov_mean(TRA,INV)(x,w,na.rm=na.rm)

#' Geometric mean
#'
#' Obtain the geometric mean of an an object. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=log` and `INV=exp`.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The geometric mean.
#' @examples 
#' geometric_mean(c(1,10,100))
#' @family Kolmogorov means
#' @export
geometric_mean <- make_Kolmogorov_mean(log, exp)

#' @rdname geometric_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_geometric_mean(c(1,10,100),c(0,1,1))
#' @export
weighted_geometric_mean <- make_weighted_Kolmogorov_mean(log, exp)

#' @rdname geometric_mean
#' @return The geometric standard deviation.
#' @examples 
#' 
#' geometric_sd(c(1,10,100))
geometric_sd <- make_Kolmogorov_mean(log, exp, sd)

#' @rdname geometric_mean
weighted_geometric_sd <- make_weighted_Kolmogorov_mean(log, exp, FUNC=sd)

#' CAGR
#' 
#' Obtain the compound annual growth rate (CAGR) of an object.
#' @param x An object consisting of a number of percentage changes over a number of years.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The CAGR.
#' @examples 
#' CAGR(c(-10,60))
#' @family Kolmogorov means
#' @export
CAGR <- make_Kolmogorov_mean(\(z) log(1+z/100), \(z) (exp(z)-1)*100)

#' @rdname CAGR
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_CAGR(c(-10,60),c(1,2))
#' @export
weighted_CAGR <- make_weighted_Kolmogorov_mean(\(z) log(1+z/100), \(z) (exp(z)-1)*100)

#' Harmonic mean
#'
#' Obtain the harmonic mean of an an object, also known as the subcontrary mean. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=\(z) 1/z` and `INV=\(z) 1/z`.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The harmonic mean.
#' @examples 
#' harmonic_mean(c(1,10,100))
#' @family Kolmogorov means
#' @export
harmonic_mean <- make_Kolmogorov_mean(\(z) 1/z, \(z) 1/z)

#' @rdname harmonic_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_harmonic_mean(c(1,10,100),c(0,1,1))
#' @export
weighted_harmonic_mean <- make_weighted_Kolmogorov_mean(\(z) 1/z, \(z) 1/z)

#' Quadratic mean
#'
#' Obtain the quadratic mean of an an object, also known as the root mean square or RMS. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=\(z) z*z` and `INV=sqrt`.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The quadratic mean.
#' @examples 
#' quadratic_mean(1:3)
#' @family Kolmogorov means
#' @export
quadratic_mean <- make_Kolmogorov_mean(\(z) z*z, \(z) sqrt(z))

#' @rdname quadratic_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples 
#' 
#' weighted_quadratic_mean(1:4,c(1,1,1,0))
#' @export
weighted_quadratic_mean <- make_weighted_Kolmogorov_mean(\(z) z*z, sqrt)

#' Cubic mean
#'
#' Obtain the cubic mean of an an object. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=\(z) z*z*z` and `INV=\(z) z^(1/3)`.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The cubic mean.
#' @examples 
#' cubic_mean(1:3)
#' @family Kolmogorov means
#' @export
cubic_mean <- make_Kolmogorov_mean(\(z) z*z*z, \(z) z^(1/3))

#' @rdname cubic_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_cubic_mean(1:4,c(1,1,1,0))
#' @export
weighted_cubic_mean <- make_weighted_Kolmogorov_mean(\(z) z*z*z, \(z) z^(1/3))

#' Square root mean
#'
#' Obtain the square root mean of an an object. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=sqrt` and `INV=\(z) z*z`. This is the opposite mean to the quadratic mean/RMS.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The square root mean.
#' @examples 
#' sqrt_mean(1:3)
#' @family Kolmogorov means
#' @export
sqrt_mean <- make_Kolmogorov_mean(sqrt, \(z) z*z)

#' @rdname sqrt_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#'  
#' weighted_sqrt_mean(1:4,c(1,1,1,0))
#' @export
weighted_sqrt_mean <- make_weighted_Kolmogorov_mean(sqrt, \(z) z*z)

#' Log semiring mean
#'
#' Obtain the log semiring mean of an an object. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=exp` and `INV=log`. This is the opposite mean to the geometric mean. The `LogSumExp` function is similar but uses `FUNC=sum` rather than `FUNC=mean`.
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The log semiring mean.
#' @examples 
#' log_semiring_mean(1:3)
#' @family Kolmogorov means
#' @export
log_semiring_mean <- make_Kolmogorov_mean(exp, log)

#' @rdname log_semiring_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples
#' 
#' weighted_log_semiring_mean(1:4,c(1,1,1,0))
#' @export
weighted_log_semiring_mean <- make_weighted_Kolmogorov_mean(exp, log)

#' @rdname log_semiring_mean
#' @return The LogSumExp value of the object.
#' @examples
#' 
#' LogSumExp(1:4)
#' @export
LogSumExp <- make_Kolmogorov_mean(exp, log, sum)

#' @rdname log_semiring_mean
#' @return The weighted LogSumExp value of the object.
#' @examples
#' 
#' weighted_LogSumExp(1:4,c(1,1,1,0))
#' @export
weighted_LogSumExp <- make_weighted_Kolmogorov_mean(exp, log, sum)

#' Holder mean functions
#'
#' Create a Holder mean function using a function factory. Also known as a power mean. This function is derived from the function factory `make_Kolmogorov_mean` using the transformation functions `FUN=\(z) z^p` and `INV=\(z) z^(1/p)` except for certain values of `p` where predefined functions are returned instead.
#' @param p The power for the transformation.
#' @details For very high or low values of `p` overflow or underflow during internal calculations may generate unexpected results. For repeat use use associated factory function and apply it, rather than the single-function-call equivalent.
#' @return A function with parameters `(x, na.rm=FALSE)` to calculate a Holder mean of an object, derived from `make_Kolmogorov_mean` except for the following cases. For `p==-Inf` the `min` function is returned. For `p==-1` the `harmonic_mean` function is returned. For `p==0` the `geometric_mean` function is returned. For `p==0.5` the `sqrt_mean` function is returned. For `p==1` the `mean` function is returned. For `p==2` the `quadratic_mean` function is returned. For `p==3` the `cubic_mean` function is returned. For `p==Inf` the `max` function is returned.
#' @examples 
#' geomean <- make_Holder_mean(0)
#' geomean(c(1,10,100))
#' @family Kolmogorov means
#' @export
make_Holder_mean <- function(p) {
  if (p==-Inf) return(min)
  else if (p==-1) return(harmonic_mean)
  else if (p==0) return(geometric_mean)
  else if (p==0.5) return(sqrt_mean)
  else if (p==1) return(mean)
  else if (p==2) return(quadratic_mean)
  else if (p==3) return(cubic_mean)
  else if (p==Inf) return(max)
  make_Kolmogorov_mean(\(z) z^p, \(z) z^(1/p))
}

#' @rdname make_Holder_mean
#' @return A function with parameters `(x, w, na.rm=FALSE)` to calculate a weighted Holder mean of an object, derived from `make_weighted_Kolmogorov_mean` except for the following cases. For `p==-Inf` the `weighted_min` function is returned. For `p==-1` the `weighted_harmonic_mean` function is returned. For `p==0` the `weighted_geometric_mean` function is returned. For `p==0.5` the `weighted_sqrt_mean` function is returned. For `p==1` the `weighted_mean` function is returned. For `p==2` the `weighted_quadratic_mean` function is returned. For `p==3` the `weighted_cubic_mean` function is returned. For `p==Inf` the `weighted_max` function is returned.
#' @examples 
#' 
#' geomeanw <- make_weighted_Holder_mean(0)
#' geomeanw(c(1,10,100,1000),c(1,1,1,0))
#' @export
make_weighted_Holder_mean <- function(p) {
  if (p==-Inf) return(weighted_min)
  else if (p==-1) return(weighted_harmonic_mean)
  else if (p==0) return(weighted_geometric_mean)
  else if (p==0.5) return(weighted_sqrt_mean)
  else if (p==1) return(weighted_mean)
  else if (p==2) return(weighted_quadratic_mean)
  else if (p==3) return(weighted_cubic_mean)
  else if (p==Inf) return(weighted_max)
  make_weighted_Kolmogorov_mean(\(z) z^p, \(z) z^(1/p))
}

#' @rdname make_Holder_mean
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return The Holder mean, calculated by a function returned by `make_Holder_mean`.
#' @examples
#' 
#' Holder_mean(c(1,10,100),0)
#' @export
Holder_mean <- function(x,p,na.rm=FALSE) make_Holder_mean(p)(x,na.rm=na.rm)

#' @rdname make_Holder_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @return The weighted Holder mean, calculated by a function returned by `make_weighted_Holder_mean`.
#' @examples 
#' 
#' weighted_Holder_mean(c(1,10,100,1000),c(1,1,1,0),0)
#' @export
weighted_Holder_mean <- function(x,w,p,na.rm=FALSE) make_weighted_Holder_mean(p)(x,w,na.rm=na.rm)
