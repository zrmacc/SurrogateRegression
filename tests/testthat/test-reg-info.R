test_that("Regression information.", {
  
  # Data.
  t <- c(1, 1, NA)
  s <- c(1, NA, 1)
  x <- c(1, 1, 1)
  sigma <- diag(c(1, 2))
  lambda <- diag(c(1, 0.5))
  
  # Partition data.
  data_part <- PartitionData(t, s, x)
  
  # Observed.
  info <- RegInfo(data_part, sigma)
  
  # Information for beta.
  ibb <- lambda[1, 1] * 1 + 1 / sigma[1, 1]
  
  # Information for alpha.
  iaa <- lambda[2, 2] * 1 + 1 / sigma[2, 2]
  
  # Cross information (note sigma is diagonal).
  iba <- 0
  
  expect_equal(as.numeric(info$Ibb), ibb)
  expect_equal(as.numeric(info$Iaa), iaa)
  expect_equal(as.numeric(info$Iba), iba)
    
})