#' Ditribution with a Power Law
#'
#' Density, distribution function, quantile function and random
#' generation for a continuous distribution with the density
#' \code{(pow+1)*(x-min)^pow/(max-min)^(pow+1)} for \code{x}
#' in the range \code{[min,max]} and \code{pow > -1}.
#' 
#' @param x vector of values where to evaluate the denisty or CDF.
#' @param p vector of probabilities.
#' @param pow degree of the power law.
#' @param n number of observations. If \code{length(n) > 1}, the length is taken to be the number required.
#' @param min minimum value of the support of the distribution.
#' @param max maximum value of the support of the distribution.
#' @returns \code{dpower} gives the density, \code{ppower} gives the cumulative
#' distribution function (CDF), \code{qpower} gives the quantile function
#' (i.e., the inverse of the CDF), and \code{rpower} generates random numbers.
#'
#' The length of the result is determined by \code{n} for \code{rpower}, and is
#' the length of \code{x} or \code{p} for the other functions.
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
