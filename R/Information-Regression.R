# Purpose: Calculate information for regression parameters.
# Updated: 2020-11-28.

#' Regression Information
#'
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param sigma Target-surrogate covariance matrix.
#' @param as_matrix Return as an information matrix? If FALSE, returns a list.
#' @return List containing the information matrix for beta (Ibb), the
#'   information matrix for alpha (Iaa), and the cross information (Iba).
RegInfo <- function(data_part, sigma, as_matrix = FALSE) {
  
  # Dimensions.
  p <- data_part$Dims$p
  q <- data_part$Dims$q
  r <- data_part$Dims$r
  n0 <- data_part$Dims$n0
  n1 <- data_part$Dims$n1
  n2 <- data_part$Dims$n2

  # Inverse covariance matrix.
  sigma_inv <- matInv(sigma)
  
  # -----------------------------------
  
  # Information for beta.
  ibb <- array(0, dim = c(q, q))
  if (n0 > 0) {
    ibb <- ibb + sigma_inv[1, 1] * data_part$IPs$X0tX0
  }
  if (n2 > 0) {
    ibb <- ibb + data_part$IPs$X2tX2 / sigma[1, 1]
  }
  
  # -----------------------------------

  # Information for alpha.
  iaa <- array(0, dim = c(r, r))
  if (n0 > 0) {
    iaa <- iaa + sigma_inv[2, 2] * data_part$IPs$Z0tZ0
  }
  if (n1 > 0) {
    iaa <- iaa + data_part$IPs$Z1tZ1 / sigma[2, 2]
  }
  
  # -----------------------------------

  # Cross information.
  iba <- array(0, dim = c(q, r))
  if (n0 > 0) {
    iba <- sigma_inv[1, 2] * data_part$IPs$X0tZ0
  }
  
  # -----------------------------------
  
  # Output.
  if (as_matrix) {
    out <- rbind(cbind(ibb, iba), cbind(t(iba), iaa))
  } else {
    out <- list(
      "Ibb" = ibb,
      "Iaa" = iaa,
      "Iba" = iba
    )
  }
  
  return(out)
}
