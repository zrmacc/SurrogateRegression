# Purpose: ECM update for covariance parameters.
# Updated: 2020-11-28.

#' Covariate Update
#'
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param b0 Previous target regression coefficient.
#' @param a0 Previous surrogate regression coefficient.
#' @param b1 Current target regression coefficient.
#' @param a1 Current surrogate regression coefficient.
#' @param sigma0 Initial target-surrogate covariance matrix.
#' @return ECM update of the target-surrogate covariance matrix.

CovUpdate <- function(data_part, b0, a0, b1, a1, sigma0) {
  
  # Dimensions.
  n <- data_part$Dims$n
  n0 <- data_part$Dims$n0
  n1 <- data_part$Dims$n1
  n2 <- data_part$Dims$n2
  
  # Structure as matrices.
  b0 <- matrix(b0, ncol = 1)
  a0 <- matrix(a0, ncol = 1)
  b1 <- matrix(b1, ncol = 1)
  a1 <- matrix(a1, ncol = 1)
  
  # Inverse covariance matrix.
  sigma_inv0 <- matInv(sigma0)

  # Residual vectors.
  et <- es <- c()
  
  # -----------------------------------
  
  # Complete cases.
  if (n0 > 0) {
    
    # Target.
    et0 <- data_part$Complete$t0 - MMP(data_part$Complete$X0, b1)
    et <- c(et, et0)
    
    # Surrogate.
    es0 <- data_part$Complete$s0 - MMP(data_part$Complete$Z0, a1)
    es <- c(es, es0)
  }
  
  # -----------------------------------
  
  # Target missing.
  if (n1 > 0) {
    
    # Surrogate.
    es1 <- data_part$TMiss$s1 - MMP(data_part$TMiss$Z1, a1)
    es <- c(es, es1)
    
    # Target.
    w1 <- (sigma0[1, 2] / sigma0[2, 2])
    et1 <- MMP(data_part$TMiss$X1, b0 - b1) + 
      w1 * (data_part$TMiss$s1 - MMP(data_part$TMiss$Z1, a0))
    et <- c(et, et1)
  }
  
  # -----------------------------------
  
  # Surrogate missing.
  if (n2 > 0) {
    
    # Target.
    et2 <- data_part$SMiss$t2 - MMP(data_part$SMiss$X2, b1)
    et <- c(et, et2)
    
    # Surrogate.
    w2 <- (sigma0[2, 1] / sigma0[1, 1])
    es2 <- MMP(data_part$SMiss$Z2, a0 - a1) + 
      w2 * (data_part$SMiss$t2 - MMP(data_part$SMiss$X2, b0))
    es <- c(es, es2)
  }

  # Residual matrix.
  E <- cbind(et, es)
  
  # Expected outer product.
  V <- matIP(E, E) + 
    n1 * diag(c(1 / sigma_inv0[1, 1], 0)) + 
    n2 * diag(c(0, 1 / sigma_inv0[2, 2]))

  ## Update sigma.
  sigma1 <- (V / n)
  return(sigma1)
}
