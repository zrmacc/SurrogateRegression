# Purpose: GLS update for regression parameters.
# Updated: 2020-11-28.

#' Regression Update
#'
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param sigma Target-surrogate covariance matrix.
#' @return List containing the generalized least squares estimates of beta and
#'   alpha.

RegUpdate <- function(data_part, sigma) {
  
  # Dimensions.
  p <- data_part$Dims$p
  q <- data_part$Dims$q
  r <- data_part$Dims$r
  n0 <- data_part$Dims$n0
  n1 <- data_part$Dims$n1
  n2 <- data_part$Dims$n2

  # Inverse target-surrogate covariance.
  sigma_inv <- matInv(sigma)

  # Information matrix for gamma = c(beta, alpha).
  reg_info <- RegInfo(data_part = data_part, sigma = sigma)
  igg <- rbind(cbind(reg_info$Ibb, reg_info$Iba), cbind(t(reg_info$Iba), reg_info$Iaa))
  iggi <- matInv(igg)

  # Calculate cross products.
  b <- array(0, dim = c(q, 1))
  a <- array(0, dim = c(r, 1))
  
  # Complete cases.
  if (n0 > 0) {
    b <- b + sigma_inv[1, 1] * data_part$IPs$X0tT0 + sigma_inv[1, 2] * data_part$IPs$X0tS0
    a <- a + sigma_inv[2, 1] * data_part$IPs$Z0tT0 + sigma_inv[2, 2] * data_part$IPs$Z0tS0
  }
  
  # Target missing.
  if (n1 > 0) {
    a <- a + data_part$IPs$Z1tS1 / sigma[2, 2]
  }
  
  # Surrogate missing.
  if (n2 > 0) {
    b <- b + data_part$IPs$X2tT2 / sigma[1, 1]
  }

  # GLS estimate of gamma.
  gls <- as.numeric(MMP(iggi, rbind(b, a)))

  # Output.
  out <- list(
    "b" = gls[1:q],
    "a" = gls[(q + 1):p]
  )
  return(out)
}
