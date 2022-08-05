# Purpose: Fitting procedure for bivariate normal regression via LS
# Updated: 2022-08-05


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
FitBNLS <- function(t, s, X, sig = 0.05) {
  
  # Partition subjects.
  data_part <- PartitionData(t, s, X)
  q <- data_part$Dims$q

  # Stage 1 regression.
  fit_1 <- fitOLS(
    y = data_part$Orig$s, 
    X = data_part$Orig$X
  )
  alpha <- fit_1$Beta
  ss_cov <- fit_1$V
  rm(fit_1)

  # Stage 2 regression.
  Z0 <- cbind(data_part$Complete$s0, data_part$Complete$X0)
  fit_2 <- fitOLS(y = data_part$Complete$t0, X = Z0)
  zeta <- fit_2$Beta
  tt_cov_inv <- fit_2$V
  rm(fit_2)

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
    method = "LS",
    b = beta,
    a = alpha,
    sigma = sigma,
    sig = sig
  )
  return(out)
}
