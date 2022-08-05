test_that("Expectation maximization Wald test.", {
  
  # Data
  withr::local_seed(101)
  n <- 1e3
  X <- cbind(1, rnorm(n))
  Z <- cbind(1, rnorm(n))
  data <- rBNR(
    X = X,
    Z = Z,
    b = c(1, 0),
    a = c(-1, 0),
    t_miss = 0.1,
    s_miss = 0.1
  )
  

  # Test 1st coefficient.
  score_test1 <- ScoreBNEM(
    t = data[, 1],
    s = data[, 2],
    X = X,
    Z = Z,
    is_zero = c(TRUE, FALSE)
  )
  expect_equal(as.numeric(score_test1["p"]), 0)

  # Test 2nd coefficient.
  score_test2 <- ScoreBNEM(
    t = data[, 1],
    s = data[, 2],
    X = X,
    Z = Z,
    is_zero = c(FALSE, TRUE)
  )
  
  expect_gt(score_test2["p"], 0.05)
  
  # Covariates passed as dataframe.
  expect_error({
    ScoreBNEM(
      t = data[, 1],
      s = data[, 2],
      X = data.frame(X),
      Z = NULL,
      is_zero = c(TRUE, FALSE)
    )
  }, NA)
  
    
})