test_that("mboot parameter", {
  set.seed(100)
  n <- 1000
  data <- runif(n)
  max.tau <- function(x) { return(x) }
  max.statistic <- function(data, indices) { return(max(data[indices])) }
  # empty data
  expect_error(mboot(c(), statistic, 0), "no data.*")
  # negative R
  expect_error(mboot(data, statistic, R = -1), ".*>0.*")
  # call to mboot which should work
  boot.out <- mboot(data, max.statistic, R = 1000, m = sqrt(n), replace = FALSE)

  # mboot.ci
  # confidence level out of bounds
  expect_error(mboot.ci(boot.out, conf = 1.01, max.tau), "conf.*")

  # user provides tau_n directly instead of a function
  expect_error(mboot.ci(boot.out, conf = 0.95, tau = max.tau(n)), ".*function")

  # mboot.ci call expected to work
  cis <- mboot.ci(boot.out, conf = 0.95, tau = max.tau)
})

test_that("estimate m", {
  set.seed(100)
  n <- 1000
  data <- runif(n)
  max.tau <- function(x) { return(x) }
  max.statistic <- function(data, indices) { return(max(data[indices])) }

  # unknown method
  expect_error(estimate.m(data, max.statistic, tau = max.tau, R = 1000, method = "unknown"),
               "unsupported.*")

  # bickel
  # choosing a bad q value
  expect_error(estimate.m(data, max.statistic,
                          tau = max.tau, R = 1000, method = "bickel", params = list(q = 0.9999)), ".*q value.*")
  set.seed(99)
  n <- 100
  data <- runif(n)
  m.bickel <- estimate.m(data, max.statistic, tau = max.tau, R = 1000, method = "bickel")
  expect_lte(m.bickel, n)
  expect_gte(m.bickel, 2, n)
  expect(m.bickel == 4, "m.bickel not 4")

  # goetze
  m.goetze <- estimate.m(data, max.statistic, tau = max.tau, R = 1000, method = "goetze")
  expect_lte(m.goetze, n)
  expect_gte(m.goetze, 2, n)

  # politis
  m.politis <- estimate.m(data, max.statistic, tau = max.tau, R = 1000, method = "politis")
  expect_lte(m.politis, n)
  expect_gte(m.politis, 2, n)

  # sherman will be tested in its own test file
})

test_that("estimate tau", {
  set.seed(100)
  # using mean
  data <- runif(10000)
  mean.statistic <- function(data, indices) { return(mean(data[indices])) }
  mean.tau <- log(estimate.tau(data, mean.statistic, R = 1000, replace = FALSE)(length(data)), length(data))
  expect_lt(mean.tau, 0.6)
  expect_gt(mean.tau, 0.4)
  mean.tau.quantile <- log(estimate.tau(data, mean.statistic, R = 1000, replace = FALSE, method = "quantile")(length(data)), length(data))
  expect_equal(mean.tau, mean.tau.quantile, tolerance = 0.05)
})
