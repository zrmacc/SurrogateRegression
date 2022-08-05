test_that("Least squares Wald test.", {
  
  # Data.
  withr::local_seed(101)
  n <- 1e3
  X <- cbind(1, stats::rnorm(n))
  data <- rBNR(
    X = X,
    Z = X,
    b = c(1, 0),
    a = c(-1, 0),
    t_miss = 0.1,
    s_miss = 0.0
  )

  # Test 1st coefficient.
  wald_test1 <- WaldBNLS(
    t = data[, 1],
    s = data[, 2],
    X = X,
    is_zero = c(TRUE, FALSE)
  )
  
  expect_equal(as.numeric(wald_test1["p"]), 0)

  # Test 2nd coefficient.
  wald_test2 <- WaldBNLS(
    t = data[, 1],
    s = data[, 2],
    X = X,
    is_zero = c(FALSE, TRUE)
  )
  expect_gt(wald_test2["p"], 0.05)
  
  # Covariates passed as dataframe.
  expect_error({
    WaldBNLS(
      t = data[, 1],
      s = data[, 2],
      X = X,
      is_zero = c(TRUE, FALSE)
    )
  }, NA)
    
})