test_that("sherman ci", {
  set.seed(100)
  max.statistic <- function(data, indices) {return(max(data[indices]))}
  data <- runif(100)
  boot.out <- mboot(data, max.statistic,m = sqrt(length(data)))
  ci <- mboot.ci(boot.out, conf = 0.95, types = "sherman")$sherman
  expect_gte(ci[2],ci[1])
  expect_lte(ci[1], ci[2])
})


test_that("sherman double bootstrap", {
  skip("Skipped for runtime reasons")
  set.seed(100)
  max.statistic <- function(data, indices) {return(max(data[indices]))}
  data <- runif(10)
  estimated.m <- estimate.m(data, max.statistic, method = "sherman")
  expect_lte(estimated.m, length(data))
  expect_gte(estimated.m, 2)
})
