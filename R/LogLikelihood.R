# Purpose: Function to calculate the observed data log likelihood.
# Updated: 2022-08-04

#' Observed Data Log Likelihood
#'
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param b Target regression coefficient.
#' @param a Surrogate regression coefficient.
#' @param sigma Target-surrogate covariance matrix.
#' @return Observed data log likelihood.
ObsLogLik <- function(data_part, b, a, sigma) {
  
  # Ensure beta and alpha are matrices.
  b <- as.matrix(b, ncol = 1)
  a <- as.matrix(a, ncol = 1)
  
  # Dimensions.
  n0 <- data_part$Dims$n0
  n1 <- data_part$Dims$n1
  n2 <- data_part$Dims$n2
  
  # -----------------------------------

  # Contribution of complete cases.
  t0 <- data_part$Complete$t0
  s0 <- data_part$Complete$s0
  X0 <- data_part$Complete$X0
  Z0 <- data_part$Complete$Z0

  resid0 <- cbind(t0 - MMP(X0, b), s0 - MMP(Z0, a))
  loglik0 <- n0 * log(matDet(sigma)) + tr(MMP(matInv(sigma), matIP(resid0, resid0)))
  
  # -----------------------------------

  # Contribution of subjects with target missingness.
  s1 <- data_part$TMiss$s1
  Z1 <- data_part$TMiss$Z1

  resid1 <- matrix(c(s1 - MMP(Z1, a)), ncol = 1)
  loglik1 <- n1 * log(sigma[2, 2]) + as.numeric(matIP(resid1, resid1)) / sigma[2, 2]
  
  # -----------------------------------

  # Contribution of subjects with surrogate missingness.
  t2 <- data_part$SMiss$t2
  X2 <- data_part$SMiss$X2

  resid2 <- matrix(c(t2 - MMP(X2, b)), ncol = 1)
  loglik2 <- n2 * log(sigma[1, 1]) + as.numeric(matIP(resid2, resid2)) / sigma[1, 1]
  
  # -----------------------------------

  # Final log likelihood
  loglik <- -0.5 * (loglik0 + loglik1 + loglik2)

  # Return
  return(loglik)
}
