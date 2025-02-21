#' m-Out-of-n Bootstrap Implementation
#'
#' Generate \code{R} bootstrap replicates of the given \code{statistic} applied to the \code{data}.
#' Sampling can be done with or without replacement.
#' The subsample size m can either be chosen directly or estimated with [estimate.m()].
#'
#' @param data The data to be bootstrapped. If it is multidimensional, each row is considered as one observation passed to the \code{statistic}.
#' @param statistic A function returing the statistic of interest. It must take two arguments. The first argument passed will be the original data, the second
#' will be a vector of indicies. Any further arguments can be passed through the \code{...} argument.
#' @param R The number of bootstrap replicates.
#' @param m The subsampling size.
#' @param replace Whether sampling should be done with replacement or without replacement (the default).
#' @param ... Additional parameters to be passed to the \code{statistic}.
#' @return The returned value is an object of the class \code{"mboot"} containing the following components:
#' \itemize{
#'  \item t0: The observed value of \code{statistic} applied to the \code{data}.
#'  \item t: A matrix with \code{R} rows where each is a bootstrap replicate of the result of calling \code{statistic}.
#'  \item m,n: Selected subsample size and data size.
#'  \item data: The \code{data} passed to \code{mboot}.
#'  \item statistic: The \code{statistic} passed to \code{mboot}.
#'  \item replace: Whether the bootstrap replicates were done with or without replacement.
#' }
#'
#' @details
#'
#' \code{m} needs to be a numeric value meeting the condition \code{2<=m<=n}.
#' It must be chosen such that m goes to infinity as n goes to infinits,
#' but the ratio m/n must go to zero.
#' The m-out-of-n Bootstrap without replacement, known as subsampling, was introduced by Politis and Romano (1994).
#' @importFrom methods hasArg
#' @export
#' @examples
#' data <- runif(1000)
#' estimate.max <- function(data, indices) {return(max(data[indices]))}
#' boot.out <- mboot(data, estimate.max, R = 1000, m = 2*sqrt(NROW(data)), replace = FALSE)
#'
#' @seealso mboot.ci estimate.m estimate.tau
#'
#' @references Politis D.N. and Romano J.P. (1994) Large sample confidence regions
#' based on subsamples under minimal assumptions. \emph{The Annals of Statistics}, 22(4):2031-2050, \doi{10.1214/aos/1176325770}
#' @keywords ~htest ~nonparametric
mboot <- function(data, statistic, m, R = 1000, replace = FALSE, ...) {
  n <- NROW(data)

  if (is.null(n) | n == 0)
    stop("no data in mboot call")
  if (R <= 0)
    stop("R needs to be >0")

  if (!hasArg(m)) {
    stop("m not provided")
  }
  if (m > n || m < 2)
    stop("m must be in [2,n]")

  t0 <- statistic(data, 1:n, ...)

  indices <- replicate(R, sample(1:n, size = m, replace = replace))
  t <- apply(indices, MARGIN = 2, FUN = function(indices) statistic(data, indices, ...))

  boot.out <- list(t0 = t0, t = t, m = m, n = n, data = data, statistic = statistic, replace = replace)
  class(boot.out) <- "mboot"
  return(boot.out)
}


#' m-Out-of-n Bootstrap Confidence Intervals
#'
#' Estimates the confidence interval using the methods provided by \code{types}.
#' \code{tau} must be a function that calculates teh scaling factor
#' tau(n) for a given n. If \code{tau} is not provided, it is estimated
#' with \code{estimate.tau} using the default settings of this function.
#'
#' @param boot.out The simulated bootstrap distribution from the \code{mboot} call.
#' @param conf The confidence level.
#' @param tau Function that returns the scaling factor tau in dependence of n. If \code{NULL}, \code{estimate.tau} is used to estimate \code{tau}.
#' @param types The types of confidence intervals to be calculated. The value can be 'all' for all types, or a
#' subset of \code{c("basic", "norm", "sherman")}.
#' @param ... When \code{tau} is omitted, the additional parameters are passed to \code{statistic} during estimation of \code{tau}.
#'
#' @returns A list of confidence intervals for the given types.
#'
#'
#' @details
#' As estimating the scaling factor tau(n) can be unreliable, it is recommended
#' to explicitly provide \code{tau}. Otherwise it is estimated with
#' \code{estimate.tau}. To specify additional arguments for
#' \code{estimate.tau}, call this function directly and use its return value
#' as \code{tau} argument. For the type \code{sherman}, \code{tau} is not
#' needed and its value is ignored.
#' 
#' The following methods to compute teh confidence intervals are supported
#' through the parameter \code{type}:
#' 
#' \describe{\item{basic:}{
#' This method works for all estimators and computes the interval directly from the quantiles of the m-out-of-n bootstrap distribution.}
#' \item{norm:}{
#' This method only works for normally distributed estimators. It estimates the variance with the m-out-of-n bootstrap and then computes te interval with the quantiles of teh standard normal distribution.}
#' \item{sherman:}{
#' This method does not scale the interval with tau(m)/tau(n) and thus is too wide. To avoid over-coverage, this is compensated by centering it randomly around the point estimators of one of the m-out-of-n bootstrap samples. Although this results on average in the nominal coverage probability, the interval is less accurate than the other intervals and should be used only as a last resort if the scaling factor tau is neither known, nor estimatable.}}
#'
#' @examples
#' data <- runif(1000)
#' estimate.max <- function(data, indices) {return(max(data[indices]))}
#' tau <- function(n){n} # convergence rate (usually sqrt(n), but n for max) 
#' boot.out <- mboot(data, estimate.max, R = 1000, m = 2*sqrt(NROW(data)), replace = FALSE)
#' cis <- mboot.ci(boot.out, 0.95, tau, c("all"))
#' ci.basic <- cis$basic
#' print(ci.basic)
#'
#' @importFrom stats quantile
#' @importFrom stats var
#' @importFrom stats ecdf
#' @importFrom stats qnorm
#'
#' @seealso mboot estimate.tau
#'
#' @references Politis D.N. and Romano J.P. (1994) Large sample confidence regions
#' based on subsamples under minimal assumptions. \emph{The Annals of Statistics}, 22(4):2031-2050, \doi{10.1214/aos/1176325770}
#' @references Sherman M. and Carlstein E. (2004) Confidence intervals based on estimators with unknown rates of convergence.
#' \emph{Computional statistics & data analysis}, 46(1):123-136.
#' @references Dalitz C. and Lögler M. (2024) moonboot: An R Package Implementing m-out-of-n Bootstrap Methods \doi{10.48550/arXiv.2412.05032}
#' @keywords ~htest
#' @export
mboot.ci <- function(boot.out, conf = 0.95, tau = NULL, types = "all", ...) {
  if (conf > 1 | conf < 0)
    stop("conf must be between 0,1")
  if (is.null(boot.out) || !inherits(boot.out, "mboot"))
    stop("you need to call mboot first")

  if ("all" %in% types)
    types <- c("basic", "norm", "sherman")

  if ("basic" %in% types | "norm" %in% types) { # tau is not needed if sherman is the only type
    if (is.null(tau))
      tau <- estimate.tau(boot.out$data, boot.out$statistic, R = 1000, replace = FALSE, ...)
    if (!is.function(tau))
      stop("tau must be a function")
  }
  alpha <- 1 - conf

  t <- boot.out$t
  t0 <- boot.out$t0

  m <- boot.out$m
  n <- boot.out$n
  cis <- list()

  if ("basic" %in% types) {
    q <- quantile((t - t0), c(1 - alpha / 2, alpha / 2)) * tau(m) / tau(n)
    upQ <- q[[1]]
    lowQ <- q[[2]]
    perc <- c(t0 - upQ, t0 - lowQ)
    cis$basic <- perc
  }
  if ("norm" %in% types) {
    z <- qnorm(1 - alpha / 2)
    sigma <- sqrt(var(t)) * tau(m) / tau(n)
    cis$norm <- c(t0 - z * sigma, t0 + z * sigma)
  }
  if ("sherman" %in% types) {
    point.est <- boot.out$t0
    centered.statistics <- boot.out$t - point.est
    G <- ecdf(centered.statistics)
    q <- quantile(G, c(1 - alpha / 2, alpha / 2))
    t.m <- sample(boot.out$t, 1) # Sherman centers the interval arround a random area of size m
    cis$sherman <- c(t.m - q[[1]], t.m - q[[2]])
  }
  return(cis)
}

#' Estimating a Subsample Size m
#'
#' Estimates \code{m} using the selected \code{method}.
#' Additional parameters can be passed to the underlying methods using \code{params}.
#' It is also possible to pass parameters to the statistic using '...'.
#'
#' @param data The data to be bootstrapped.
#' @param statistic The estimator of the parameter.
#' @param tau The convergence rate.
#' @param R The amount of bootstrap replicates. Must be a positive integer.
#' @param method The method to be used, one of \code{c("goetze","bickel","politis", "sherman")}.
#' @param replace If the sampling should be done with replacement. Setting this value to true requires a sufficient smooth estimator.
#' @param min.m Minimum subsample size to be tried. Should be the minimum size for which the statistic make sense.
#' @param params Additional parameters to be passed to the internal functions, see details for more information.
#' @param ... Additional parameters to be passed to the statistic.
#'
#' @return Subsampling size \code{m} choosen by the selected method.
#'
#' @details
#' The different methods have different parameters. Therefore, this wrapper method has been given the \code{params} parameter, which can be used to
#' pass method-specific arguments to the underlying methods. The specific parameters are described below.
#' Most of the provided methods need \code{tau}. If not provided, it will be estimated using
#' \code{estimate.tau}. Note that method 'sherman' is using an alternative approach without using the scalation factor and
#' therefore \code{tau} will not be computed if selecting 'sherman' as method. Any non \code{NULL} values will be ignored when
#' selecting the method 'sherman'.
#'
#' Possible methods are:
#' 
#' \describe{\item{goetze:}{
#' The method from Goetze and Rackauskas is based on minimizing the distance between the
#' CDF of the bootstrap distributions of different subsampling sizes 'm'.
#' As distance measurement the 'Kolmogorov distance' is used.
#' The method uses the pairs 'm' and 'm/2' to be minimized.
#' As this would involve trying out all combinations of 'm' and 'm/2' this method has a running time of order Rn^2.
#' To reduce the runtime in practical use, \code{params} can be used to pass a \code{goetze.interval}, which is a
#' list of the smallest and largest value for m to try.}
#' \item{bickel:}{
#' This method works similary to the previous one. The difference here is that the subsample sizes to be
#' compared are consecutive subsample sizes generated by \code{q^j*n} for \code{j = seq(2,n)} and a chosen \code{q} value between
#' zero and one.
#' The parameter \code{q} can be selected using \code{params}. The default value is \code{q=0.75}, as suggested in the corresponding paper.}
#' \item{politis:}{
#' This method is also known as the 'minimum volatility method'. It is based on the idea that there
#' should be some range for subsampling sizes, where its choice has little effect on the estimated confidence points.
#' The algorithm starts by smoothing the endpoints of the intervals and then calculates the standard deviation.
#' The \code{h.ci} parameter is used to select the number of neighbors used for smoothing.
#' The \code{h.sigma} parameter is the number of neighbors used in the standard deviation calculation.
#' Both parameters can be set by using \code{params}.
#' Note that the \code{h.*} neigbors from each side are used.
#' To use five elements for smoothing, \code{h.ci} should therefore be set to 2.}
#' \item{sherman:}{
#' This method is based on a 'double-bootstrap' approach.
#' It tries to estimate the coverage error of different subsampling sizes and chooses the subsampling
#' size with the lowest one.
#' As estimating the coverage error is highly computationally intensive, it is not practical to try all m values.
#' Therefore, the \code{beta} parameter can be used to control which \code{m} values are tried. The values
#' are then calculated by \code{ms = n^beta}. The default value is a sequence between 0.3 and 0.9 out of 15 values.
#' This parameter can be set using \code{params}.}}
#'
#' @examples
#' data <- runif(1000)
#' estimate.max <- function(data, indices) {return(max(data[indices]))}
#' tau <- function(n){n} # convergence rate (usually sqrt(n), but n for max) 
#' choosen.m <- estimate.m(data, estimate.max, tau, R = 1000, method = "bickel")
#' print(choosen.m)
#'
#'
#' @seealso mboot estimate.tau
#'
#' @references Götze F. and Rackauskas A. (2001) Adaptive choice of bootstrap sample sizes.
#' \emph{Lecture Notes-Monograph Series}, 36(State of the Art in Probability and Statistics):286-309
#' @references Bickel P.J. and Sakov A. (2008) On the choice of m in the m out of n bootstrap and confidence bounds for extrema.
#' \emph{Statistic Sinica}, 18(3):967-985.
#' @references Politis D.N. et al. (1999)
#' \emph{Subsampling}, Springer, New York.
#' @references Sherman M. and Carlstein E. (2004) Confidence intervals based on estimators with unknown rates of convergence.
#' \emph{Computional statistics & data analysis}, 46(1):123-136.
#'
#' @importFrom stats filter
#' @importFrom stats sd
#' @export

estimate.m <- function(data, statistic, tau = NULL, R = 1000, replace = FALSE, min.m = 3, method = "bickel", params = NULL, ...) {
  if (is.null(tau) & method != "sherman") {
    tau <- estimate.tau(data, statistic, R = 1000, replace = FALSE, ...)
  }

  if ("goetze" == method) {
    return(estimate.m.goetze(data, statistic, tau, R, replace, min.m, params, ...))
  }
  else if ("bickel" == method) {
    return(estimate.m.bickel(data, statistic, tau, R, replace, min.m, params, ...))
  }else if ("politis" == method) {
    return(estimate.m.politis(data, statistic, tau, R, replace, min.m, params, ...))
  }else if ("sherman" == method) {
    return(estimate.m.sherman(data, statistic, R, replace, min.m, params, ...))
  }else {
    stop("unsupported method to estimate m")
  }
}

estimate.m.goetze <- function(data, statistic, tau, R, replace = FALSE, min.m, params, ...) {
  n <- NROW(data)
  params.default.values <- list(goetze.interval = c(4, n))

  if (!hasArg(params) || is.null(params)) {
    params <- params.default.values
  } else if (!is.list(params)) {
    stop("params needs to be a list")
  }
  # warning if params has unknown parameters which are not in params.default.values
  unexpected.params <- setdiff(names(params), names(params.default.values))
  if (length(unexpected.params) > 0) {
    warning(paste("estimate.m.goetze: The following parameters are unknown and will be ignored: ", paste(unexpected.params, collapse = ", ")))
  }
  # min, max value of search interval
  goetze.interval <- if ("goetze.interval" %in% names(params)) params[["goetze.interval"]] else params.default.values[["goetze.interval"]]

  # calculate t0 once
  t0 <- statistic(data, 1:n, ...)

  # Calculate the distance between m and m/2
  calcDistance <- function(bout.m, bout.mh) {
    m <- bout.m$m
    m.half <- bout.mh$m

    m.t <- bout.m$t
    mh.t <- bout.mh$t

    m.t <- tau(m) * (m.t - t0)
    mh.t <- tau(m.half) * (mh.t - t0)
    F <- ecdf(m.t)
    F2 <- ecdf(mh.t)
    allres <- c(m.t, mh.t)
    distance <- max(abs(F(allres) - F2(allres)))
    return(distance)
  }

  # using an (optional) search interval to speed up the process
  a <- max(2 * min.m, goetze.interval[1] + goetze.interval[1] %% 2) # limiting search interval to min.m
  b <- min(n, goetze.interval[2])
  # filtering m values which are not in the search interval
  m.values <- seq(a, b, by = 2)
  m.values <- m.values[m.values >= 2 * min.m] # it is needed to calculate m/2, therefore m needs to be at least 2*min.m
  m.values <- m.values[m.values <= n]
  required.m.values <- unique(c(m.values, ceiling(m.values / 2)))
  required.m.values <- sort(required.m.values)
  required.boots <- lapply(required.m.values, function(m)
    mboot(data, statistic, R, m = m, replace = replace, ...))
  # iterate over m.values and calculate distances between corresponding required.boots m and m/2.
  distances <- sapply(m.values, function(m)
    calcDistance(required.boots[[which(required.m.values == m)]], required.boots[[which(required.m.values == ceiling(m / 2))]]))

  m <- m.values[which.min(distances)]
  return(m)
}

estimate.m.bickel <- function(data, statistic, tau, R, replace = FALSE, min.m, params, ...) {
  n <- NROW(data)
  params.default.values <- c(q = 0.75)

  if (!hasArg(params) || is.null(params)) {
    params <- params.default.values
  } else if (!is.list(params)) {
    stop("params needs to be a list")
  }
  # warning if params has unknown parameters which are not in params.default.values
  unexpected.params <- setdiff(names(params), names(params.default.values))
  if (length(unexpected.params) > 0) {
    warning(paste("estimate.m.bickel: The following parameters are unknown and will be ignored: ", paste(unexpected.params, collapse = ", ")))
  }
  q <- if ("q" %in% names(params)) params[["q"]] else params.default.values["q"]

  calcDistance <- function(m, mp) {
    m.boot <- mboot(data, statistic, R, m = m, replace = replace, ...)
    mp.boot <- mboot(data, statistic, R, m = mp, replace = replace, ...)

    m.t <- m.boot$t
    m.star.t <- mp.boot$t
    m.t <- tau(m) * (m.t - m.boot$t0)
    m.star.t <- tau(mp) * (m.star.t - m.boot$t0)
    F <- ecdf(m.t)
    F2 <- ecdf(m.star.t)
    allres <- c(m.t, m.star.t)
    distance <- max(abs(F(allres) - F2(allres)))
    return(distance)
  }

  j <- seq(2, n)
  mj <- q^j * n
  mj <- ceiling(mj)

  mj <- mj[!duplicated(ceiling(n / mj))]
  mj <- mj[mj >= min.m & mj <= n]

  if (length(mj) < 3) {
    stop("There were not enough m values calculated. Please try a different q value.")
  }

  distances <- sapply(1:length(mj[-length(mj)]), function(m) calcDistance(mj[m], mj[m + 1]))
  m.distances <- cbind(ceiling(mj[seq_along(distances)]), distances)
  # get the right minimum of distances
  min.index <- length(distances) - which.min(rev(distances)) + 1
  return(m.distances[min.index, 1])
}

estimate.m.politis <- function(data, statistic, tau, R, replace = FALSE, min.m, params, ...) {

  params.default.values <- list(h.sigma = 3, h.ci = 3, conf = 0.95)

  if (!hasArg(params) || is.null(params)) {
    params <- params.default.values
  } else if (!is.list(params)) {
    stop("params needs to be a list")
  }
  # warning if params has unknown parameters which are not in params.default.values
  unexpected.params <- setdiff(names(params), names(params.default.values))
  if (length(unexpected.params) > 0) {
    warning(paste("estimate.m.politis: The following parameters are unknown and will be ignored: ", paste(unexpected.params, collapse = ", ")))
  }

  h.sigma <- if ("h.sigma" %in% names(params)) params[["h.sigma"]] else params.default.values[["h.sigma"]]
  h.ci <- if ("h.ci" %in% names(params)) params[["h.ci"]] else params.default.values[["h.ci"]]
  conf <- if ("conf" %in% names(params)) params[["conf"]] else params.default.values[["conf"]]

  n <- NROW(data)
  m.values <- seq(min.m, n * 0.5, length.out = 50) # limit to n/2
  m.values <- round(m.values)
  m.values <- unique(m.values)
  bouts <- lapply(m.values, function(m) mboot(data, statistic, R, m = m, replace = replace, ...))
  cis <- lapply(bouts, function(bout) mboot.ci(bout, conf = conf, tau = tau, c("basic"))$basic)

  # 2. as approximation was used, smooth lower and upper endpoints seperately
  # using a running mean
  rmean.size <- h.ci * 2 + 1 # both sides
  cis.upper <- lapply(cis, function(x) x[1])
  cis.lower <- lapply(cis, function(x) x[2])
  smoothed.cis.upper <- filter(cis.upper, sides = 2, rep(1 / rmean.size, rmean.size))
  smoothed.cis.lower <- filter(cis.lower, sides = 2, rep(1 / rmean.size, rmean.size))
  # edge cases will have NA as value - calculate them manually
  last.na <- (length(smoothed.cis.lower) - h.ci + 1):length(smoothed.cis.lower)
  # no mirroring is used as this could negativly impact the variance calculation
  for (i in 1:h.ci) {
    # first.na values going from 1:h.ci, using i directly
    smoothed.cis.lower[i] <- mean(unlist(cis.lower[1:(h.ci + i)]))
    smoothed.cis.upper[i] <- mean(unlist(cis.upper[1:(h.ci + i)]))

    smoothed.cis.upper[last.na[i]] <- mean(unlist(cis.upper[(last.na[i] - h.ci):length(smoothed.cis.lower)]))
    smoothed.cis.lower[last.na[i]] <- mean(unlist(cis.lower[(last.na[i] - h.ci):length(smoothed.cis.lower)]))
  }

  VI.upper <- c()
  VI.lower <- c()
  VI <- c()

  # 3. for each m compute a volatility Index as standard deviation in their neighborhood of m
  for (i in 1:length(smoothed.cis.upper)) {
    # check if edge case
    if (i - h.sigma < 1 || i + h.sigma > length(smoothed.cis.upper)) {
      VI[i] <- Inf # ignoring edge cases
      next
    }

    # sd of neighborhood m-h.sigma to m+h.sigma
    VI.upper[i] <- sd(smoothed.cis.upper[max(1, i - h.sigma):min(length(smoothed.cis.upper), i + h.sigma)])
    VI.lower[i] <- sd(smoothed.cis.lower[max(1, i - h.sigma):min(length(smoothed.cis.lower), i + h.sigma)])
    VI[i] <- VI.upper[i] + VI.lower[i]
  }
  # Pick the value m corresponding to the smallest VI as m.star
  m.star <- m.values[which.min(VI)]

  return(m.star)
}


#' Estimating the convergence rate
#'
#' This function estimates the convergence rate of the bootstrap estimator
#' and returns it as a function of the form \code{tau_n = n^a}, where \code{n} is the input parameter.
#'
#' @param data The data to be bootstrapped.
#' @param statistic The estimator of the parameter.
#' @param R Amount of bootstrap replicates used to estimate tau.
#' @param replace If sampling should be done with replacement.
#' @param min.m Minimal subsampling size used to estimate tau. Should be set to the minimum size for which the statistic makes sense.
#' @param beta The tested subsample sizes m are \code{n^beta}.
#' @param method Method to estimate tau, can be one of \code{c("variance", "quantile")}.
#' @param ... Additional parameters to be passed to the \code{mboot} function.
#'
#' @return A function for the square root of the convergence rate of the variance, i.e., \code{f(n) = tau_n}. This function can directly be passed to \code{mboot.ci}.
#'
#' @details There are two methods to choose from, \code{variance} and \code{quantile}.
#' The provided \code{beta} values are used to select subsample sizes \code{m} by using \code{ms = n^beta}.
#' Note that the choice of the \code{beta} values can impact the accuracy of the estimated \code{tau} (Dalitz & Lögler, 2024).
#' For each selected subsample size \code{m} a bootstrap with \code{R} replications is performed.
#' The method 'variance' then fits a linear function to log(variance) of the bootstrap statistics as function of log(m).
#' The method 'quantile' averages over multiple quantile ranges Q and fits a linear function to log(Q) as a function of log(m).
#'
#'
#' @examples
#' data <- runif(1000)
#' estimate.max <- function(data, indices) {return(max(data[indices]))}
#' estimated.tau <- estimate.tau(data, estimate.max)
#' boot.out <- mboot(data, estimate.max, R = 1000, m = 2*sqrt(NROW(data)), replace = FALSE)
#' cis <- mboot.ci(boot.out, 0.95, estimated.tau, c("all"))
#' ci.basic <- cis$basic
#' print(ci.basic)
#'
#' @seealso mboot.ci
#'
#' @references Bertail P. et al. (1999) On subsampling estimators with unknown rate of convergence.
#' \emph{Journal of the American Statistical Association}, 94(446):568-579.
#' @references Politis D.N. et al. (1999)
#' \emph{Subsampling}, Springer, New York.
#' @references Dalitz, C, and Lögler, F. (2024)
#' \emph{moonboot: An R Package Implementing m-out-of-n Bootstrap Methods}.
#' \doi{10.48550/arXiv.2412.05032}
#'
#' @importFrom stats lm
#' @importFrom stats coef
#' @export
estimate.tau <- function(data, statistic, R = 1000, replace = FALSE, min.m = 3, beta = seq(0.2, 0.7, length.out = 5), method = "variance", ...) {
  if ("variance" == method) {
    return(estimate.tau.variance(data, statistic, R = R, replace = replace, min.m = min.m, beta = beta, ...))
  }else if ("quantile" == method) {
    return(estimate.tau.quantile(data, statistic, R = R, replace = replace, min.m = min.m, bs = beta, ...))
  }
}

# Estimates tau using the variance method suggested by Bertail et al. (1999)
estimate.tau.variance <- function(data, statistic, R = 1000, replace = FALSE, min.m, beta = seq(0.2, 0.7, length.out = 5), ...) {
  n <- NROW(data)
  beta <- pmin(1, pmax(0, beta)) # force them to be in bounds
  # b from the paper is called m for consistency with the rest of the code
  m <- n^beta
  m <- m[m >= min.m]
  m.boots <- lapply(m, function(m) mboot(data, statistic, R = R, m = ceiling(m), replace = replace, ...))
  V.m <- lapply(m.boots, function(boot) var(boot$t))
  V.m <- unlist(V.m)
  y <- -log(V.m) / 2
  x <- log(m)
  a <- lm(y ~ x)
  a <- coef(a)
  return(function(n) n^a[[2]])
}

# Estimate tau using the method suggested by Polits et al. (1999).
estimate.tau.quantile <- function(data, statistic, R = 1000, replace = FALSE, min.m, bs, ...) {
  n <- NROW(data)
  J <- 5
  I <- 5
  j <- seq(1, J, by = 1)
  t <- numeric(2 * J)

  seq.one <- seq(from = 0.75, by = 0.05, length.out = J)
  seq.two <- seq(from = 0.25, by = -0.05, length.out = J)
  for (i in j) {
    t[2 * i] <- seq.one[i]
    t[2 * i - 1] <- seq.two[i]
  }
  if (!hasArg(bs))
    bs <- seq(0.2, 0.7, length.out = I)
  ms <- floor(n^bs)
  ms <- unique(ms)
  ms <- ms[ms >= min.m & ms <= n]
  I <- length(ms)
  boot.outs <- lapply(ms, function(m) mboot(data, statistic, R = R, m = m, replace = replace, ...))

  get.log.difference <- function(i, j) {
    # tau is set to 1 according to the formula
    mout <- boot.outs[[i]]
    ci <- quantile(mout$t, c(t[2 * j], t[2 * j - 1]), na.rm = TRUE)
    first.lower <- ci[1]
    second.lower <- ci[2]
    difference <- first.lower - second.lower
    log.dif <- log(difference)
    return(log.dif)
  }

  # multi dimensional matrix usage: y.ij[i,j]
  y.ij <- matrix(0, nrow = I, ncol = J)
  for (i in 1:I) {
    for (j in 1:J) {
      y.ij[i, j] <- get.log.difference(i, j)
      if (is.infinite(y.ij[i, j])) {
        y.ij[i, j] <- 0
      }
    }
  }
  yi <- numeric(I)
  for (i in 1:I) {
    yi[i] <- 1 / J * sum(y.ij[i,])
  }
  y.bar <- mean(yi)
  log.bar <- 1 / I * sum(log(ms))
  beta.denom <- (log(ms) - log.bar)^2

  beta.numerator <- (yi - y.bar) * (log(ms) - log.bar)

  beta.IJ <- -sum(beta.numerator) / sum(beta.denom)
  return(function(x) x^beta.IJ)
}

