#' Mean of the Shorthest Half
#'
#' Calculates the mean of the data points in the shortest interval containing half of the data.
#' This function is ready to be used as a statistic in [mboot()] function.
#'
#' @param data the data as a vector.
#' @param indices the selected indices of \code{data}, by default \code{seq_along(data)}.
#' @returns The mean of the data points in the shortest interval containing half of the data.
#'
#' @examples
#' data <- rnorm(100)
#' shorth(data)
#' shorth(data, sample(1:100, size = 100))
#'
#' # Calculating a CI for shorth using [mboot()]
#' data <- rnorm(100)
#' boot.out <- mboot(data, shorth, m = sqrt(length(data)))
#' basic.ci <- mboot.ci(boot.out, conf =0.95, tau = function(n) return(n^(1/3)), types = "basic")$basic
#'
#' @references Andrews D.F. et al. (1972) \emph{Robust Estimates of Location} Princeton University Press, Princeton.
#' @export
shorth <- function(data, indices = NULL) {
  if (is.null(indices))
    indices <- seq_along(data)
  x <- data[indices]
  x <- sort(x)
  n <- length(x)

  mid <- ceiling(n / 2)
  # generate all possible intervals
  if (n %% 2 == 0)
    left.indices <- 1:(mid + 1)
  else
    left.indices <- 1:(mid)
  right.indices <- mid:n
  interval.lengths <- x[right.indices] - x[left.indices]
  min.dist <- which.min(interval.lengths)
  l <- left.indices[min.dist]
  r <- right.indices[min.dist]
  return(mean(x[l:r]))
}