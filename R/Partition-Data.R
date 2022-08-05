# Purpose: Function to partition subjects by outcome missingness pattern.
# Updated: 2020-10-25.


#' Partition Data by Outcome Missingness Pattern.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @return List containing these components:
#' \itemize{
#'   \item `Orig` original data.
#'   \item `Dims` dimensions and names.
#'   \item `Complete`, data for complete cases.
#'   \item `TMiss`, data for subjects with target missingness.
#'   \item `SMiss`, data for subjects with surrogate missingness.
#'   \item `IPs`, inner products. 
#' }
#' @export
#' @examples
#' # Generate data.
#' n <- 1e3
#' X <- rnorm(n)
#' Z <- rnorm(n)
#' data <- rBNR(X = X, Z = Z, b = 1, a = -1)
#' data_part <- PartitionData(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z
#' )
PartitionData <- function(t, s, X, Z = NULL) {
  
  # Ensure input structures are matrices.
  t <- matrix(t, ncol = 1)
  s <- matrix(s, ncol = 1)
  if (!is.matrix(X)) {X <- as.matrix(X)}
  
  # If no surrogate model matrix is provided, 
  # adopt the target model matrix.
  if (is.null(Z)) {
    Z <- X
  } else {
    if (!is.matrix(Z)) {Z <- as.matrix(Z)}
  }
  
  # Output structure.
  out <- list()
  
  # -------------------------------------------------------
  
  out$Orig <- list()
  out$Orig$t <- t
  out$Orig$s <- s
  out$Orig$X <- X
  out$Orig$Z <- Z
  
  # -------------------------------------------------------
  
  # Design matrix dimensions.
  q <- ncol(X)
  r <- ncol(Z)
  p <- q + r

  out$Dims <- list()
  out$Dims$p <- p
  out$Dims$q <- q
  out$Dims$r <- r
  
  # Names.
  x_names <- colnames(X)
  if (is.null(x_names)) {
    x_names <- paste0("x", seq_len(ncol(X)))
  }
  z_names <- colnames(Z)
  if (is.null(z_names)) {
    z_names <- paste0("z", seq_len(ncol(Z)))
  }
  
  out$Dims$x_names <- x_names
  out$Dims$z_names <- z_names
  
  # -------------------------------------------------------

  # Observation indicators.
  n <- nrow(t)
  R <- array(data = 0, dim = c(n, 2))
  R[, 1] <- 1 * (!is.na(t))
  R[, 2] <- 1 * (!is.na(s))

  out$Dims$n <- n
  out$Inds <- list()
  out$Inds$R <- R
  
  # -------------------------------------------------------

  # Partition complete cases.
  is_case_0 <- (R[, 1] == 1) & (R[, 2] == 1)
  t0 <- t[is_case_0, , drop = FALSE]
  s0 <- s[is_case_0, , drop = FALSE]
  X0 <- X[is_case_0, , drop = FALSE]
  Z0 <- Z[is_case_0, , drop = FALSE]
  n0 <- length(t0)

  # Storage.
  out$Complete <- list()
  out$Complete$t0 <- t0
  out$Complete$s0 <- s0
  out$Complete$X0 <- X0
  out$Complete$Z0 <- Z0
  out$Dims$n0 <- n0
  
  # -------------------------------------------------------

  # Target missing, surrogate observed.
  is_case_1 <- (R[, 1] == 0) & (R[, 2] == 1)
  t1 <- t[is_case_1, , drop = FALSE]
  s1 <- s[is_case_1, , drop = FALSE]
  X1 <- X[is_case_1, , drop = FALSE]
  Z1 <- Z[is_case_1, , drop = FALSE]
  n1 <- length(t1)

  # Storage.
  out$TMiss <- list()
  out$TMiss$t1 <- t1
  out$TMiss$s1 <- s1
  out$TMiss$X1 <- X1
  out$TMiss$Z1 <- Z1
  out$Dims$n1 <- n1
  
  # -------------------------------------------------------

  # Surrogate missing, surrogate observed.
  is_case_2 <- (R[, 1] == 1) & (R[, 2] == 0)
  t2 <- t[is_case_2, , drop = FALSE]
  s2 <- s[is_case_2, , drop = FALSE]
  X2 <- X[is_case_2, , drop = FALSE]
  Z2 <- Z[is_case_2, , drop = FALSE]
  n2 <- length(t2)

  # Storage.
  out$SMiss <- list()
  out$SMiss$t2 <- t2
  out$SMiss$s2 <- s2
  out$SMiss$X2 <- X2
  out$SMiss$Z2 <- Z2
  out$Dims$n2 <- n2
  
  # -------------------------------------------------------

  # Inner products.
  out$IPs <- list()

  # Complete cases.
  out$IPs$X0tX0 <- matIP(X0, X0)
  out$IPs$X0tZ0 <- matIP(X0, Z0)
  out$IPs$Z0tZ0 <- matIP(Z0, Z0)
  out$IPs$X0tT0 <- matIP(X0, t0)
  out$IPs$X0tS0 <- matIP(X0, s0)
  out$IPs$Z0tT0 <- matIP(Z0, t0)
  out$IPs$Z0tS0 <- matIP(Z0, s0)

  # Target missing.
  out$IPs$Z1tZ1 <- matIP(Z1, Z1)
  out$IPs$Z1tS1 <- matIP(Z1, s1)

  # Surrogate missing.
  out$IPs$X2tX2 <- matIP(X2, X2)
  out$IPs$X2tT2 <- matIP(X2, t2)

  # Output.
  return(out)
}
