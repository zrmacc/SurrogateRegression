# Purpose: Calculate information matrix for covariance parameters.
# Updated: 2020-11-28

#' Covariance Information Matrix
#'
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param sigma Target-surrogate covariance matrix.
#' @return 3x3 Numeric information matrix for the target variance,
#'   target-surrogate covariance, and surrogate variance.

CovInfo <- function(data_part, sigma) {
  
  # Dimensions.
  n0 <- data_part$Dims$n0
  n1 <- data_part$Dims$n1
  n2 <- data_part$Dims$n2

  # Lambda.
  sigma_inv <- matInv(sigma)

  out <- array(0, dim = c(3, 3))
  
  # Complete cases.
  if (n0 > 0) {
    obs_info <- array(0, dim = c(3, 3))
    obs_info[1, 1] <- sigma_inv[1, 1]^2
    obs_info[2, 2] <- 2 * (sigma_inv[1, 2]^2 + sigma_inv[1, 1] * sigma_inv[2, 2])
    obs_info[3, 3] <- sigma_inv[2, 2]^2
    obs_info[1, 2] <- obs_info[2, 1] <- 2 * sigma_inv[1, 1] * sigma_inv[1, 2]
    obs_info[2, 3] <- obs_info[3, 2] <- 2 * sigma_inv[1, 2] * sigma_inv[2, 2]
    obs_info[1, 3] <- obs_info[3, 1] <- sigma_inv[1, 2]^2
    out <- out + 0.5 * n0 * obs_info
  }
  
  # Target missing, only surrogate observed.
  if (n1 > 0) {
    out[3, 3] <- out[3, 3] + 0.5 * n1 / (sigma[2, 2]^2)
  }
  
  # Surrogate missing, only target observed.
  if (n2 > 0) {
    out[1, 1] <- out[1, 1] + 0.5 * n2 / (sigma[1, 1]^2)
  }
  return(out)
}
