knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 6,
  fig.height = 4
)
library(SurrogateRegression)

set.seed(100)
n <- 800L
X <- cbind(1, rnorm(n))   # target design (intercept + 1 covariate)
Z <- cbind(1, rnorm(n))   # surrogate design
Y <- rBNR(
  X = X,
  Z = Z,
  b = c(1, 1),            # target coefficients
  a = c(-1, -1),          # surrogate coefficients
  t_miss = 0.2,           # 20% target missing
  s_miss = 0.1            # 10% surrogate missing
)
head(Y, 10)

fit <- FitBNR(
  t = Y[, "Target"],
  s = Y[, "Surrogate"],
  X = X,
  Z = Z
)

show(fit)

coef(fit)
coef(fit, type = "Target")

head(residuals(fit), 5)
head(residuals(fit, type = "Target"), 5)

# Residual covariance of (target, surrogate)
vcov(fit, type = "Outcome")
# Information matrix for regression coefficients (inverse = asymptotic cov of coefs)
info_reg <- vcov(fit, type = "Regression")
dim(info_reg)

# Test that the first target coefficient (intercept) is zero
test_intercept <- TestBNR(
  t = Y[, "Target"],
  s = Y[, "Surrogate"],
  X = X,
  Z = Z,
  is_zero = c(TRUE, FALSE),
  test = "Wald"
)
test_intercept

# Test that the second target coefficient is zero
test_slope <- TestBNR(
  t = Y[, "Target"],
  s = Y[, "Surrogate"],
  X = X,
  Z = Z,
  is_zero = c(FALSE, TRUE),
  test = "Wald"
)
test_slope

part <- PartitionData(
  t = Y[, "Target"],
  s = Y[, "Surrogate"],
  X = X,
  Z = Z
)
names(part)
part$Dims$n0   # complete cases
part$Dims$n1   # target missing, surrogate observed
part$Dims$n2   # surrogate missing, target observed
