#' Ditribution with a Power Law
#'
#' Density, distribution function, quantile function and random
#' generation for a continuous distribution with the density
#' (pow+1)*(x-min)^pow for x in `[min,max]` and pow > -1.
#' 
#' @param x vector of values where to evaluate the denisty or CDF.
#' @param p vector of probabilities.
#' @param pow degree of the power law.
#' @param n number of observations. If 'length(n) > 1', the length is taken to be the number reqired.
#' @param min minimum value of the support of the distribution.
#' @param max maximum value of the support of the distribution.
#' @returns ‘dpower’ gives the density, ‘ppower’ gives the cumulative
#' distribution function (CDF), ‘qpower’ gives the quantile function
#' (i.e., the inverse of the CDF), and ‘rpower’ generates random numbers.
#'
#' The length of the result is determined by ‘n’ for ‘rpower’, and is
#' the length of x or p for the other functions.
#' @name distPower
NULL

#' @rdname distPower
#' @export
dpower <- function(x, pow, min=0, max=1) {
  stopifnot(pow > -1)
  ifelse(x>min & x<max, (pow+1)*(x-min)^pow/(max-min)^(pow+1), 0)
}

#' @rdname distPower
#' @export
ppower <- function(x, pow, min=0, max=1) {
  stopifnot(pow > -1)
  ifelse(x<min, 0, ifelse(x>max, 1, (x-min)^(pow+1)/(max-min)^(pow+1)))
}

#' @rdname distPower
#' @export
qpower <- function(p, pow, min=0, max=1) {
  stopifnot(pow > -1)
  stopifnot(p>=0 & p<=1)
  min + p^(1/(pow+1)) * (max-min)
}

#' @rdname distPower
#' @importFrom stats runif
#' @export
rpower <- function(n, pow, min=0, max=1) {
  stopifnot(pow > -1)
  if(length(n) > 1)
    n <- length(n)
  x <- runif(n)
  min + x^(1/(pow+1)) * (max-min)
}
