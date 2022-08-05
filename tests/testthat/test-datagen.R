test_that("Data generation.", {
  
  n <- 4
  x <- rep(1, n)
  z <- x
  b <- 1
  a <- 1
  
  # Calculation of design matrix.
  obs <- rBNR(x, z, b, a, include_residual = FALSE)
  exp <- array(1.0, dim = c(n, 2))
  expect_equal(obs, exp, ignore_attr = TRUE)
  
  # Target missingness only.
  y <- rBNR(x, z, b, a, t_miss = 0.25, include_residuals = FALSE)
  t <- y[, 1]
  s <- y[, 2]
  obs <- c(sum(is.na(t)), sum(is.na(s)))
  exp <- c(1, 0)
  expect_equal(obs, exp)
  
  # Surrogate missingness only.
  y <- rBNR(x, z, b, a, s_miss = 0.50, include_residuals = FALSE)
  t <- y[, 1]
  s <- y[, 2]
  obs <- c(sum(is.na(t)), sum(is.na(s)))
  exp <- c(0, 2)
  expect_equal(obs, exp)
  
  # Bilateral missingness: 1.
  y <- rBNR(
    x, z, b, a, t_miss = 0.25, s_miss = 0.25, include_residuals = FALSE)
  t <- y[, 1]
  s <- y[, 2]
  obs <- c(sum(is.na(t)), sum(is.na(s)), sum(is.na(t) & is.na(s)))
  exp <- c(1, 1, 0)
  expect_equal(obs, exp)
  
  # Bilateral missingness: 2.
  y <- rBNR(
    x, z, b, a, t_miss = 0.50, s_miss = 0.50, include_residuals = FALSE)
  t <- y[, 1]
  s <- y[, 2]
  obs <- c(sum(is.na(t)), sum(is.na(s)), sum(is.na(t) & is.na(s)))
  exp <- c(2, 2, 0)
  expect_equal(obs, exp)
  
  # Bilateral missingness: 3.
  expect_error({
    rBNR(
      x, z, b, a, t_miss = 0.75, s_miss = 0.50, include_residuals = FALSE)
  })
  
})