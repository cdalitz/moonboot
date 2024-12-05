# internal functions suggested by Sherman and Carlstein

# double bootstrap approach to choose a m value
estimate.m.sherman <- function(data, statistic, R = 1000, replace = FALSE, min.m, params, ...) {
  params.default.values <- list(beta = seq(0.3, 0.9, length.out = 15), conf = 0.95)

  if (!hasArg(params) | is.null(params)) {
    params <- params.default.values
  } else if (!is.list(params)) {
    stop("params needs to be a list")
  }
  # warning if params has unknown parameters which are not in params.default.values
  unexpected.params <- setdiff(names(params), names(params.default.values))
  if (length(unexpected.params) > 0) {
    warning(paste("estimate.m.sherman: The following parameters are unknown and will be ignored: ", paste(unexpected.params, collapse = ", ")))
  }
  # min, max value of search interval
  beta <- if ("beta" %in% names(params)) params[["beta"]] else params.default.values[["beta"]]
  conf <- if ("conf" %in% names(params)) params[["conf"]] else params.default.values[["conf"]]

  threshold <- 0.05
  n <- NROW(data)
  ms <- round(n^beta)
  ms <- unique(ms)
  ms <- ms[ms > min.m & ms < n]
  model.point.est <- statistic(data, 1:n)
  cov.ps <- lapply(ms, function(m) sherman.estimate.covp(data, statistic, m, model.point.est, conf, R, replace = replace,...))
  cov.ps <- as.numeric(unlist(cov.ps))

  # choose beta where two elements in a row are smaller than choosen threshold
  ms.index <- 1
  for (ms.i in 1:(length(ms) - 3)) {
    if (cov.ps[ms.i + 1] >= threshold | cov.ps[ms.i + 2] >= threshold) {
      ms.index <- ms.i
    }
  }
  choosen.m <- round(ms[ms.index])
  return(choosen.m)
}


# internal method to estimate the coverage probability for a choosen m value
sherman.estimate.covp <- function(data, statistic, m, t0, conf, N, replace = FALSE, ...) {
  n <- NROW(data)
  generator.replace <- TRUE

  # generator function to generate samples from the data
  generator <- function(data, n) {
    if (is.vector(data)) {
      return(data[sample(1:n, size = n, replace = generator.replace)])
    } else if (is.matrix(data)) {
      return(data[sample(1:n, size = n, replace = generator.replace), , drop = FALSE])
    } else {
      stop("Unsupported data type")
    }
  }

  coverage.sum <- 0
  for (i in 1:N) {
    gen.data <- generator(data, n)
    bout <- mboot(gen.data, statistic, m = m, R = 1000, replace = replace, ...)
    ci <- mboot.ci(bout, types = "sherman", conf = conf)$sherman
    coverage.sum <- coverage.sum + ifelse(t0 >= ci[1] & t0 <= ci[2], 1, 0)
  }
  return(coverage.sum / N)
}
