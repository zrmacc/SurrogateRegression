test_that("Log likelihood calculation.", {

  # Data.
  x <- c(1, 1, 1)
  t <- c(NA, 1, 1)
  s <- c(1, NA, 1)
  data_part <- PartitionData(t, s, x, x)
  
  # Case 1.
  a <- 0
  b <- 0
  sigma <- diag(2)
  
  obs <- ObsLogLik(data_part, b = b, a = a, sigma = sigma)
  exp <- -0.5 * (2 + 1 + 1)
  expect_equal(obs, exp)
  
  # Case 2.
  a <- 1
  b <- 0
  sigma <- diag(2)
  
  obs <- ObsLogLik(data_part, b = b, a = a, sigma = sigma)
  exp <- -0.5 * (1 + 1)
  expect_equal(obs, exp)
  
  # Case 3.
  a <- 0
  b <- 0
  sigma <- 2 * diag(2)
  
  obs <- ObsLogLik(data_part, b = b, a = a, sigma = sigma)
  exp <- -0.5 * (log(4) + 2 * log(2) + (1 + 0.5 + 0.5))
  expect_equal(obs, exp)
  
})