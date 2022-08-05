test_that("Fit bivariate model via expectation maximization.", {
  
  # Data.
  withr::local_seed(101)
  n <- 1e3
  x <- rep(1, n)
  b <- 1
  a <- -1
  data <- rBNR(X = x, Z = x, b = b, a = a, t_miss = 0.1, s_miss = 0.1)
  
  # Observed.
  obs <- FitBNEM(t = data[, 1], s = data[, 2], X = x, Z = x, report = FALSE)
  obs <- obs@Regression.tab
  
  expect_equal(obs$Point[1], b, tolerance = 0.05)
  expect_equal(obs$Point[2], a, tolerance = 0.05)
  
})