# Purpose: Wald test for bivariate normal regression via EM.
# Updated: 19/08/02

#' Wald Test via Expectation Maximization.
#'
#' Performs a Wald test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @param is_zero Logical vector, with as many entries as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param init Optional list of initial parameters for fitting the null model,
#'   with one or more of the components: a0, b0, S0.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#'
#' @return A numeric vector containing the Wald statistic, the degrees of
#'   freedom, and a p-value.
#'   
#' @examples 
#' \donttest{
#' # Generate data.
#' set.seed(100)
#' n <- 1e3
#' X <- cbind(1, rnorm(n))
#' Z <- cbind(1, rnorm(n))
#' data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.1, s_miss = 0.1)
#' 
#' # Test 1st coefficient.
#' wald_test1 <- Wald.BNEM(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(TRUE, FALSE)
#' )
#' 
#' # Test 2nd coefficient.
#' wald_test2 <- Wald.BNEM(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(FALSE, TRUE)
#' )
#' }

Wald.BNEM <- function(
  t, 
  s, 
  X, 
  Z, 
  is_zero, 
  init = NULL, 
  maxit = 100, 
  eps = 1e-8, 
  report = FALSE
) {

  # Input checks.
  CheckInit(init = init)
  CheckTestSpec(is_zero = is_zero, p = ncol(X))

  # For Wald test, fit the full model.
  fit <- Fit.BNEM(
    t = t, 
    s = s, 
    X = X, 
    Z = Z, 
    b0 = init$b0, 
    a0 = init$a0, 
    sigma0 = init$sigma0, 
    maxit = maxit, 
    eps = eps, 
    report = report
  )

  # Extract information.
  reg_info <- vcov(fit, type = "Regression", inv = FALSE)
  
  # Surrogate covariates.
  q <- ncol(Z)

  # Partition regression information into the components zero under the null
  # ('key0'), and the components unconstrained under the null ('key1').
  key0 <- c(is_zero, rep(FALSE, q))
  key1 <- c(!is_zero, rep(TRUE, q))

  # Partition information.
  ibb <- reg_info[key0, key0, drop = FALSE]
  iaa <- reg_info[key1, key1, drop = FALSE]
  iba <- reg_info[key0, key1, drop = FALSE]
  
  # Efficient information.
  inv_var <- SchurC(Ibb = ibb, Iaa = iaa, Iba = iba)

  # Coefficients of interest.
  beta_hat <- coef(fit, type = "Target")$Point[is_zero]
  beta_hat <- matrix(beta_hat, ncol = 1)
  
  # Wald statistic.
  stat_wald <- as.numeric(matQF(X = beta_hat, A = inv_var))
  
  # P value.
  df <- sum(is_zero)
  pval <- pchisq(q = stat_wald, df = df, lower.tail = FALSE)
  
  # Output.
  out <- c(
    "Wald" = stat_wald, 
    "df" = df, 
    "p" = pval
  )
  return(out)
}
