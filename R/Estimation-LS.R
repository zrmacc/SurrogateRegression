# Purpose: Fitting procedure for bivariate normal regression via LS
# Updated: 2020-11-28.

#' Fit Bivariate Normal Regression Model via Least Squares
#'
#' Estimation procedure for bivariate normal regression models in which
#' only the target outcome is subject to missingness.
#'
#' The model matrix is expected in numeric format. Include an intercept if
#' required. Expand factors and interactions in advance.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Model matrix.
#' @param sig Type I error level.
#' @return An object of class 'bnr' with slots containing the estimated
#'   regression coefficients, the target-surrogate covariance matrix, the
#'   information matrices for the regression and covariance parameters, and the
#'   residuals.
#'
#' @importFrom stats coef pnorm qnorm resid
#' @export
#' 
#' @examples 
#' \donttest{
#' set.seed(100)
#' n <- 1e3
#' X <- rnorm(n)
#' data <- rBNR(X = X, Z = X, b = 1, a = -1, t_miss = 0.1, s_miss = 0.0)
#' t <- data[, 1]
#' s <- data[, 2]
#' 
#' # Model fit.
#' fit_bnls <- Fit.BNLS(
#'   t = t,
#'   s = s,
#'   X = X,
#'   sig = 0.05
#' )
#' }

Fit.BNLS <- function(t, s, X, sig = 0.05) {
  
  # Partition subjects.
  data_part <- PartitionData(t, s, X)
  q <- data_part$Dims$q

  # Stage 1 regression.
  fit.1 <- fitOLS(
    y = data_part$Orig$s, 
    X = data_part$Orig$X
  )
  alpha <- fit.1$Beta
  ss_cov <- fit.1$V
  rm(fit.1)

  # Stage 2 regression.
  Z0 <- cbind(data_part$Complete$s0, data_part$Complete$X0)
  fit.2 <- fitOLS(y = data_part$Complete$t0, X = Z0)
  zeta <- fit.2$Beta
  tt_cov_inv <- fit.2$V
  rm(fit.2)

  # Recover original parameters.
  delta <- zeta[1]
  gamma <- zeta[2:(q + 1)]
  beta <- delta * alpha + gamma
  ts_cov <- delta * ss_cov
  tt_cov <- tt_cov_inv + ss_cov * (delta)^2

  # Final covariance.
  sigma <- matrix(c(tt_cov, ts_cov, ts_cov, ss_cov), nrow = 2)
  colnames(sigma) <- rownames(sigma) <- c("Target", "Surrogate")
  
  # -----------------------------------
  
  # Output.
  out <- FormatOutput(
    data_part = data_part,
    b = beta,
    a = alpha,
    sigma = sigma,
    sig = sig
  )
  return(out)
}
