# Purpose: Wald test for bivariate normal regression via least squares.
# Updated: 2022-08-05


#' Wald Test via Least Squares.
#'
#' Performs a Wald test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Model matrix.
#' @param is_zero Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @return A numeric vector containing the Wald statistic, the degrees of
#'   freedom, and a p-value.
WaldBNLS <- function(t, s, X, is_zero) {
  
  # Input checks.
  p <- ncol(X)
  CheckTestSpec(is_zero = is_zero, p = p)

  # Model fit.
  fit <- FitBNLS(t = t, s = s, X = X)

  # Extract information.
  reg_info <- stats::vcov(fit, type = "Regression", inv = FALSE)

  # Partition regression information into the components zero under the null
  # ('key0'), and the components unconstrained under the null ('key1').
  key0 <- c(is_zero, rep(FALSE, p))
  key1 <- c(!is_zero, rep(TRUE, p))

  # Partition information.
  ibb <- reg_info[key0, key0, drop = FALSE]
  iaa <- reg_info[key1, key1, drop = FALSE]
  iba <- reg_info[key0, key1, drop = FALSE]
  
  # Efficient information.
  inv_var <- SchurC(Ibb = ibb, Iaa = iaa, Iba = iba)

  # Coefficients of interest.
  beta_hat <- stats::coef(fit, type = "Target")$Point[is_zero]
  beta_hat <- matrix(beta_hat, ncol = 1)
  
  # Statistic.
  stat_wald <- as.numeric(matQF(X = beta_hat, A = inv_var))
  
  # P value.
  df <- sum(is_zero)
  pval <- stats::pchisq(q = stat_wald, df = df, lower.tail = FALSE)
  
  # Output.
  out <- c(
    "Wald" = stat_wald, 
    "df" = df, 
    "p" = pval
  )
  return(out)
}
