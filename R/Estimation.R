# Purpose: Master estimation function for bivariate normal regression
# Updated: 2020-11-28.

#' Fit Bivariate Normal Regression Model.
#'
#' Estimation procedure for bivariate normal regression models. The EM algorithm
#' is applied if \code{s} contains missing values, or if \code{X} differs from
#' \code{Z}. Otherwise, an accelerated least squares procedure is applied.
#'
#' The target and surrogate model matrices are expected in numeric format.
#' Include an intercept if required. Expand factors and interactions in advance.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix. Defaults to \code{X}.
#' @param sig Significance level.
#' @param ... Additional arguments accepted if fitting via EM. See
#'   \code{\link{Fit.BNEM}}.
#' @return An object of class 'mnr' with slots containing the estimated regression
#'  coefficients, the target-surrogate covariance matrix, the information matrices
#'  for regression parameters, and the residuals.
#'  
#' @importFrom stats coef pnorm qnorm resid
#' @export
#' 
#' @examples
#' 
#' # Case 1: No surrogate missingness.
#' set.seed(100)
#' n <- 1e3
#' X <- rnorm(n)
#' data <- rBNR(X = X, Z = X, b = 1, a = -1, t_miss = 0.1, s_miss = 0.0)
#' t <- data[, 1]
#' s <- data[, 2]
#' 
#' # Model fit.
#' fit_bnls <- Fit.BNR(
#'   t = t,
#'   s = s,
#'   X = X
#' )
#' 
#' # Case 2: Target and surrogate missingness.
#' set.seed(100)
#' n <- 1e3
#' X <- rnorm(n)
#' Z <- rnorm(n)
#' data <- rBNR(X = X, Z = Z, b = 1, a = -1, t_miss = 0.1, s_miss = 0.1)
#' 
#' # Log likelihood.
#' fit_bnem <- Fit.BNR(
#'   t = data[, 1],
#'   s = data[, 2],
#'   X = X,
#'   Z = Z
#' )

Fit.BNR <- function(t, s, X, Z = NULL, sig = 0.05, ...) {
  
  # Determine if s contains missing values, or if Z differs from X.
  apply_em <- (sum(is.na(s)) > 0) | ((!is.null(Z)) & (!identical(X, Z)))
  
  # If so, apply EM algorithm.
  if (apply_em) {
    out <- Fit.BNEM(t = t, s = s, X = X, Z = Z, sig = sig, ...)
  } else {
    
    # Otherwise, apply the least squares procedure.
    out <- Fit.BNLS(t = t, s = s, X = X)
  }
 
   # Output
  return(out)
}
