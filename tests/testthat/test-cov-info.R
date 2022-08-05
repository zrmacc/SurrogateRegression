test_that("Covariance update.", {
  
  # Data.
  t <- c(1, 1, 1, NA)
  s <- c(1, 1, NA, 1)
  x <- c(1, 1, 1, 1)
  
  stt <- 1
  sss <- 2
  sigma <- diag(c(stt, sss))
  
  # Partition data.
  data_part <- PartitionData(t, s, x)
  
  # Observed.
  obs <- CovInfo(data_part, sigma)
  
  # Expected.
  itt <- 0.5 * (2 * (1/stt)^2 + 1 / (stt^2))
  iss <- 0.5 * (2 * (1/sss)^2 + 1 / (sss^2))
  its <- 0.5 * (2 * 2 * (1/stt) * (1/sss))
  
  expect_equal(itt, obs[1, 1])
  expect_equal(its, obs[2, 2])
  expect_equal(iss, obs[3, 3])
})