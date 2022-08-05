test_that("Covariance update.", {
  
  # Data.
  t <- c(1, 1, 1, NA)
  s <- c(1, 1, NA, 1)
  x <- c(1, 1, 1, 1)
  
  # Partition data.
  data_part <- PartitionData(t, s, x)
  
  # Residual matrix.
  e <- rbind(
    c(1, 1),
    c(1, 1),
    c(0, 1),
    c(1, 0)
  )
  exp <- t(e) %*% e + diag(c(1, 0)) + diag(c(0, 1))
  exp <- exp / nrow(e)
  obs <- CovUpdate(data_part, 0, 0, 0, 0, diag(2))
  expect_equal(obs, exp)
    
})