#' @rdname PinkAverage-package
#' @title PinkAverage
#' @description
#' Calculate various different types of averages and other related summary statistics. Individual functions for various types of means or averages are included, and function factories for creating families of means such as  Gini, Holder, and Lehmer means, both weighted and unweighted.
#' 
#' @section Author:
#' * Darren Green <darren.green@stir.ac.uk>
#' 
#' @section Function factories:
#' For most families of means, four functions are included. Function factories are preceded by `make`, e.g.
#' 
#' `f <- make_Susan_mean(params)`
#' 
#' `f(x)`
#' 
#' and
#' 
#' `f <- make_weighted_Susan_mean(params)`
#' 
#' `f(x,w)`.
#' 
#' For one-off use, equivalent functions to both make and then use the function factories are provided as
#'
#' `Susan_mean(x, params)` and
#' 
#' `weighted_Susan_mean(x, w, params)`.
#' 
#' For repeat use use associated factory function in advance and apply it, rather than the single-function-call equivalent, as this is more efficient. Most factory-produced functions support the argument `na.rm` to remove `NA` and `NaN` values before computation proceeds.
#' @section Families:
#' Several families of means are implemented, including
#' * The very general Kolmogorov mean family, including...
#' * The power mean or Holder mean family.
#' * The Lehmer mean family.
#' * The Gini mean family, which genereralises the Holder and Lehmer families.
#' * A family of truncated (including Winsorized) means.
#' * A family of n-point L estimators including midsummary functions.
#' * A family of Frechet (Karcher) means.
#' * A family of circular and modular means.
#' * A family of Bayesian averages.
#' @section Bivariate means:
#' Several families of means are implemented which operate on pairs of numbers rather than more general vectors of numbers, as follows:
#' 
#' `f <- make_Susan_mean(params)`
#' 
#' `f(x,y)`.
#' 
#' * A factory-function family for iterated means.
#' * A family of Seiffert-like means.
#' * A family of Stolarsky means, including the identric.
#' * A family of Heinz means.  
#' @section Miscellaneous means:
#' Called as follows:
#' 
#' `Susan_mean(x,params)`
#' 
#' * The mode.
#' * The pseudomedian.
#' @section Miscellaneous bivariate means:
#' Called as follows:
#' 
#' `Susan_mean(x,y)`
#' 
#' * The Heronian mean.
#' @docType package
#' @name PinkAverage
NULL
