#' @include PinkAverage-package.r

#' @rdname make_Gini_mean
#' @title Gini mean functions
#' @description Create a Gini mean function using a function factory. Gini means are a two-parameter family of means with parameters `r` and `s` which generalize the Lehmer and Holder mean families, interpolating between various means.
#' @param r,s The parameters of the mean.
#' @details For very high or low values of `p` overflow or underflow during internal calculations may generate unexpected results. For repeat use use associated factory function and apply it, rather than the single-function-call equivalent. Removal of `NA` and `NaN` values is performed by `mean`, `min`, `max` or `sum`, therefore will remove any such values introduced by manipulating the data.
#' @return A function with parameters `(x, na.rm=FALSE)` to calculate a Gini mean of an object, except for the following cases. For `s==0` the `make_Holder_mean(r)` function is returned. For `r==1` the `make_Lehmer_mean(s+1)` function is returned. For `r==-1` the `make_Lehmer_mean(s)` function is returned. For `r==0` a special-case function is returned.
#' @examples 
#' contramean <- make_Gini_mean(1,1)
#' contramean(c(1,2,3))
#' @family Gini means
#' @references 
#' Farnsworth, D and Orr, R (1986). Gini means. The American Mathematical Monthly 93: 603-607.
#' @export
make_Gini_mean <- function(r,s) {
  if (s==0) return(make_Holder_mean(r))
  else if (r==1) return(make_Lehmer_mean(s+1))
  else if (r==-1) return(make_Lehmer_mean(s))
  else if (r==0) return(function(x,na.rm=FALSE) exp(sum(x^s*log(x),na.rm=na.rm)/sum(x^s,na.rm=na.rm)))
  function(x,na.rm=FALSE) (sum(x^(r+s),na.rm=na.rm)/sum(x^s,na.rm=na.rm))^(1/r)
}

#' @rdname make_Gini_mean
#' @return A function with parameters `(x, w, na.rm=FALSE)` to calculate a weighted Gini mean of an object, with associated weights, except for the following cases. For `s==0` the `make_weighted_Holder_mean(r)` function is returned. For `r==1` the `make_weighted_Lehmer_mean(s+1)` function is returned. For `r==-1` the `make_weighted_Lehmer_mean(s)` function is returned. For `r==0` a special-case function is returned.
#' @examples 
#' 
#' contrameanw <- make_weighted_Gini_mean(1,1)
#' contrameanw(c(1,2,3,4),c(1,1,1,0))
#' @export
make_weighted_Gini_mean <- function(r,s) {
  if (s==0) return(make_weighted_Holder_mean(r))
  else if (r==1) return(make_weighted_Lehmer_mean(s+1))
  else if (r==-1) return(make_weighted_Lehmer_mean(s))
  else if (r==0) return(function(x,w,na.rm=FALSE) exp(sum(w*x^s*log(x),na.rm=na.rm)/sum(w*x^s,na.rm=na.rm)))
  function(x,w,na.rm=FALSE) (sum(w*(x^(r+s)),na.rm=na.rm)/sum(w*(x^s),na.rm=na.rm))^(1/r)   
}

#' @rdname make_Gini_mean
#' @param x An `R` object.
#' @param na.rm A logical evaluating to `TRUE` or `FALSE` indicating whether `NA` values should be stripped before the computation proceeds.
#' @return A Gini mean
#' @examples 
#' 
#' Gini_mean(c(1,2,3,4),1,1)
#' @export
Gini_mean <- function(x,r,s,na.rm=FALSE) make_Gini_mean(r,s)(x,na.rm=na.rm)

#' @rdname make_Gini_mean
#' @param w An `R` object of the same length as `x` representing a series of weights.
#' @examples 
#' weighted_Gini_mean(c(1,2,3,4,5),c(1,1,1,1,0),1,1)
#' @export
weighted_Gini_mean <- function(x,w,r,s,na.rm=FALSE) make_weighted_Gini_mean(r,s)(x,w,na.rm=na.rm)
