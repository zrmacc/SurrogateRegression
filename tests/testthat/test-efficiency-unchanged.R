# Tests to ensure efficiency optimizations do not change package behavior.
# Run with same seeds and compare outputs for determinism and optional-arg consistency.
#
# Use all.equal() + fail() for determinism checks to avoid loading waldo/diffobj,
# which can cause namespaceExport errors in some environments.

test_that("rBNR is deterministic and linearPredictors matches MMP.", {
  n <- 200
  set.seed(301)
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  b <- c(1, -0.5)
  a <- c(-1, 0.5)
  # Same seed => same data
  set.seed(301)
  y1 <- rBNR(X = X, Z = Z, b = b, a = a, t_miss = 0.1, s_miss = 0.05)
  set.seed(301)
  y2 <- rBNR(X = X, Z = Z, b = b, a = a, t_miss = 0.1, s_miss = 0.05)
  msg <- all.equal(y1, y2, tolerance = 1e-9, check.attributes = FALSE)
  if (!isTRUE(msg)) {
    testthat::fail(paste("rBNR outputs differed:", msg))
  }
  # linearPredictors matches two MMP calls
  eta_batch <- linearPredictors(X, Z, matrix(b, ncol = 1), matrix(a, ncol = 1))
  eta_sep <- cbind(MMP(X, matrix(b, ncol = 1)), MMP(Z, matrix(a, ncol = 1)))
  msg <- all.equal(eta_batch, eta_sep, tolerance = 1e-9, check.attributes = FALSE)
  if (!isTRUE(msg)) {
    testthat::fail(paste("linearPredictors vs MMP:", msg))
  }
  expect_true(TRUE)
})

test_that("FitBNR is deterministic.", {
  n <- 300
  set.seed(302)
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  set.seed(302)
  data <- rBNR(X = X, Z = Z, b = c(1, 0.5), a = c(-1, -0.5), t_miss = 0.15, s_miss = 0.1)
  fit1 <- FitBNR(t = data[, 1], s = data[, 2], X = X, Z = Z, report = FALSE)
  set.seed(302)
  data2 <- rBNR(X = X, Z = Z, b = c(1, 0.5), a = c(-1, -0.5), t_miss = 0.15, s_miss = 0.1)
  fit2 <- FitBNR(t = data2[, 1], s = data2[, 2], X = X, Z = Z, report = FALSE)
  msg <- all.equal(coef(fit1)$Point, coef(fit2)$Point, tolerance = 1e-9)
  if (!isTRUE(msg)) {
    testthat::fail(paste("coef differed:", msg))
  }
  msg <- all.equal(
    vcov(fit1, type = "Regression"),
    vcov(fit2, type = "Regression"),
    tolerance = 1e-8
  )
  if (!isTRUE(msg)) {
    testthat::fail(paste("vcov differed:", msg))
  }
  msg <- all.equal(residuals(fit1), residuals(fit2), tolerance = 1e-9)
  if (!isTRUE(msg)) {
    testthat::fail(paste("residuals differed:", msg))
  }
  expect_true(TRUE)
})

test_that("RegInfo with sigma_inv matches RegInfo without.", {
  withr::local_seed(303)
  n <- 150
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.1, s_miss = 0.1)
  part <- PartitionData(data[, 1], data[, 2], X, Z)
  sigma <- matrix(c(1.2, 0.4, 0.4, 1.1), 2, 2)
  info_null <- SurrogateRegression:::RegInfo(part, sigma, sigma_inv = NULL, as_matrix = TRUE)
  info_opt <- SurrogateRegression:::RegInfo(part, sigma, sigma_inv = matInv(sigma), as_matrix = TRUE)
  expect_equal(info_null, info_opt, tolerance = 1e-9)
})

test_that("CovInfo with sigma_inv matches CovInfo without.", {
  withr::local_seed(304)
  n <- 150
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.1, s_miss = 0.1)
  part <- PartitionData(data[, 1], data[, 2], X, Z)
  sigma <- matrix(c(1.2, 0.4, 0.4, 1.1), 2, 2)
  info_null <- SurrogateRegression:::CovInfo(part, sigma, sigma_inv = NULL)
  info_opt <- SurrogateRegression:::CovInfo(part, sigma, sigma_inv = matInv(sigma))
  expect_equal(info_null, info_opt, tolerance = 1e-9)
})

test_that("CovUpdate is deterministic and returns 2x2.", {
  withr::local_seed(305)
  n <- 100
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.2, s_miss = 0.15)
  part <- PartitionData(data[, 1], data[, 2], X, Z)
  b0 <- c(0.9, 0.1)
  a0 <- c(-0.95, 0.05)
  b1 <- c(1.02, -0.02)
  a1 <- c(-1.01, 0.01)
  sigma0 <- matrix(c(1.1, 0.35, 0.35, 1.05), 2, 2)
  out1 <- SurrogateRegression:::CovUpdate(part, b0, a0, b1, a1, sigma0)
  out2 <- SurrogateRegression:::CovUpdate(part, b0, a0, b1, a1, sigma0)
  expect_equal(out1, out2, tolerance = 1e-9)
  expect_equal(dim(out1), c(2L, 2L))
  expect_true(is.finite(out1[1, 1]) && is.finite(out1[2, 2]))
})

test_that("ObsLogLik is deterministic and matches two-call consistency.", {
  withr::local_seed(306)
  n <- 100
  X <- cbind(1, stats::rnorm(n))
  Z <- cbind(1, stats::rnorm(n))
  data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.1, s_miss = 0.1)
  part <- PartitionData(data[, 1], data[, 2], X, Z)
  b <- c(1, 0)
  a <- c(-1, 0)
  sigma <- matrix(c(1.0, 0.3, 0.3, 1.0), 2, 2)
  ll1 <- SurrogateRegression:::ObsLogLik(part, b = b, a = a, sigma = sigma)
  ll2 <- SurrogateRegression:::ObsLogLik(part, b = b, a = a, sigma = sigma)
  expect_equal(ll1, ll2, tolerance = 1e-9)
  expect_true(is.finite(ll1) && length(ll1) == 1L)
})

test_that("matInv: small matrix inverse is consistent.", {
  # 2x2 and 3x3 use inv(); larger use pinv(); both should agree for invertible matrices
  A2 <- matrix(c(1.5, 0.2, 0.2, 1.1), 2, 2)
  expect_equal(solve(A2), matInv(A2), tolerance = 1e-9, ignore_attr = TRUE)
  A3 <- matrix(c(2, 0.1, 0, 0.1, 1.5, 0.1, 0, 0.1, 1.2), 3, 3)
  expect_equal(solve(A3), matInv(A3), tolerance = 1e-8, ignore_attr = TRUE)
})
