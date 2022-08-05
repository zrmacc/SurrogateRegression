test_that("Check estimation procedure.", {
  
  # Data.
  withr::local_seed(101)
  n <- 1e3
  x <- stats::rnorm(n)
  z <- stats::rnorm(n)
  b <- 1
  a <- -1
  data <- rBNR(
    X = x,
    Z = z,
    b = b,
    a = a,
    t_miss = 0.0,
    s_miss = 0.0
  )
  
  # No missing data, different covariate sets.
  fit <- FitBNR(
    t = data[, 1],
    s = data[, 2],
    X = x,
    Z = z,
    report = FALSE
  )
  tab <- stats::coef(fit)
  expect_equal(tab$Point[1], b, tolerance = 0.05)
  expect_equal(tab$Point[2], a, tolerance = 0.05)
  expect_equal(fit@Method, "EM")  # Because X = Z, EM is required.
  
  # No missing data, same covariate set.
  fit <- FitBNR(
    t = data[, 1],
    s = data[, 2],
    X = x,
    report = FALSE
  )
  tab <- stats::coef(fit)
  expect_equal(tab$Point[1], b, tolerance = 0.05)
  expect_equal(tab$Point[2], 0, tolerance = 0.05)
  expect_equal(fit@Method, "LS")
  
  # Target missingness only.
  t_miss <- data[, 1]
  t_miss[1:5] <- NA
  fit <- FitBNR(
    t = t_miss,
    s = data[, 2],
    X = x,
    report = FALSE
  )
  tab <- stats::coef(fit)
  expect_equal(tab$Point[1], b, tolerance = 0.05)
  expect_equal(tab$Point[2], 0, tolerance = 0.05)
  expect_equal(fit@Method, "LS")
  
  # Surrogate missingness only.
  s_miss <- data[, 2]
  s_miss[1:5] <- NA
  fit <- FitBNR(
    t = data[, 1],
    s = s_miss,
    X = x,
    Z = z,
    report = FALSE
  )
  tab <- stats::coef(fit)
  expect_equal(tab$Point[1], b, tolerance = 0.05)
  expect_equal(tab$Point[2], a, tolerance = 0.05)
  expect_equal(fit@Method, "EM")
  
  # Bilateral missingness.
  fit <- FitBNR(
    t = t_miss,
    s = s_miss,
    X = x,
    Z = z,
    report = FALSE
  )
  tab <- stats::coef(fit)
  expect_equal(tab$Point[1], b, tolerance = 0.05)
  expect_equal(tab$Point[2], a, tolerance = 0.05)
  expect_equal(fit@Method, "EM")
  
  # Covariate missingness.
  x_miss <- x
  x_miss[1] <- NA
  expect_error({
    FitBNR(
      t = data[, 1],
      s = data[, 2],
      X = x_miss,
      Z = z,
      report = FALSE
    )
  })
  
  # Covariates passed as dataframes.
  fit <- FitBNR(
    t = t_miss,
    s = data[, 2],
    X = data.frame(x),
    report = FALSE
  )
  expect_equal(fit@Method, "LS")
  
  fit <- FitBNR(
    t = t_miss,
    s = data[, 2],
    X = data.frame(x),
    Z = NULL,
    report = FALSE
  )
  expect_equal(fit@Method, "LS")
  
  fit <- FitBNR(
    t = t_miss,
    s = data[, 2],
    X = data.frame(x),
    Z = data.frame(z),
    report = FALSE
  )
  expect_equal(fit@Method, "EM")
  
})