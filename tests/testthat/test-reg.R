test_that("Regression update.", {
  
  # Data.
  t <- c(1, 1, NA)
  s <- c(1, NA, 1)
  x <- c(1, 1, 1)
  sigma <- diag(c(2, 1))
  lambda <- diag(c(1/2, 1))
  
  # Partition data.
  data_part <- PartitionData(t, s, x)
  
  # Observed.
  obs <- RegUpdate(data_part, sigma)
  
  # Expected.
  igg <- diag(c(1, 2))
  g <- c(
    lambda[1, 1] * 1 + 1 / sigma[1, 1],
    lambda[2, 2] * 1 + 1 / sigma[2, 2]
  )
  exp <- solve(igg, g)
  
  expect_equal(obs$b, exp[1])
  expect_equal(obs$a, exp[2])
    
})