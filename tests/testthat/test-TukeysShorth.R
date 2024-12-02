# testing shorth on prepared data
test_that("shorth test data", {
  set.seed(100)
  data.ordered <- c(-103,-102,-101, 1, 100,101,102,103)
  data <- sample(data.ordered, replace = F) # shuffle
  shorth.mean <- shorth(data)
  expect_equal(shorth.mean, mean(data.ordered[5:8]))

  # odd number of elements
  data.ordered <- c(-103,-102,-101, 100,101,102,103)
  data <- sample(data.ordered, replace = F) # shuffle
  shorth.mean <- shorth(data)
  expect_equal(shorth.mean, mean(data.ordered[4:7]))
})


test_that("shorth normal distribution", {
  set.seed(100)
  n <- 1000
  data <- rnorm(n)
  shorth.mean <- shorth(data)

  expect_gte(shorth.mean, -0.75)
  expect_lte(shorth.mean, 0.25) 

  n <- n-1 # odd number of elements
  data <- rnorm(n)
  shorth.mean <- shorth(data)
  expect_gte(shorth.mean, -0.75)
  expect_lte(shorth.mean, 0.75)
})

