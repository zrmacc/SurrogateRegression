test_that("Test inference procedure.", {
  
  # Data
  withr::local_seed(101)
  n <- 1e3
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  data <- rBNR(
    X = X,
    Z = Z,
    b = c(1, 0),
    a = c(-1, 0),
    t_miss = 0.1,
    s_miss = 0.1
  )

  # Test 1st coefficient.
  wald_test1 <- TestBNR(
    t = data[, 1],
    s = data[, 2],
    X = X,
    Z = Z,
    is_zero = c(TRUE, FALSE),
    test = "Wald"
  )
  expect_equal(wald_test1["p"], 0, ignore_attr = TRUE)

  score_test1 <- TestBNR(
    t = data[, 1],
    s = data[, 2],
    X = X,
    Z = Z,
    is_zero = c(TRUE, FALSE),
    test = "Score"
  )
  expect_equal(score_test1["p"], 0, ignore_attr = TRUE)

  # Test 2nd coefficient.
  wald_test2 <- TestBNR(
    t = data[, 1],
    s = data[, 2],
    X = X,
    Z = Z,
    is_zero = c(FALSE, TRUE),
    test = "Wald"
  )
  expect_gt(wald_test2["p"], 0.05)

  score_test2 <- TestBNR(
    t = data[, 1],
    s = data[, 2],
    X = X,
    Z = Z,
    is_zero = c(FALSE, TRUE),
    test = "Score"
  )
  expect_gt(score_test2["p"], 0.05)
  
})