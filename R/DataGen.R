# Purpose: Data generation for bivariate normal regression model.
# Updated: 2022-08-05

# -----------------------------------------------------------------------------
# Data Generation
# -----------------------------------------------------------------------------

#' Simulate Bivariate Normal Data with Missingness
#'
#' Function to simulate from a bivariate normal regression model with outcomes
#' missing completely at random.
#'
#' @param X Target design matrix.
#' @param Z Surrogate design matrix.
#' @param b Target regression coefficient.
#' @param a Surrogate regression coefficient.
#' @param t_miss Target missingness in [0,1].
#' @param s_miss Surrogate missingness in [0,1].
#' @param sigma 2x2 target-surrogate covariance matrix.
#' @param include_residuals Include the residual? Default: TRUE.
#' @return Numeric Nx2 matrix. The first column contains the target
#'   outcome, the second contains the surrogate outcome.
#' @export
#' @examples
#' set.seed(100)
#' # Observations.
#' n <- 1e3
#' # Target design.
#' X <- cbind(1, matrix(rnorm(3 * n), nrow = n))
#' # Surrogate design.
#' Z <- cbind(1, matrix(rnorm(3 * n), nrow = n))
#' # Target coefficient.
#' b <- c(-1, 0.1, -0.1, 0.1)
#' # Surrogate coefficient.
#' a <- c(1, -0.1, 0.1, -0.1)
#' # Covariance structure.
#' sigma <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
#' # Data generation, target and surrogate subject to 10% missingness.
#' y <- rBNR(X, Z, b, a, t_miss = 0.1, s_miss = 0.1, sigma = sigma)

rBNR <- function(
  X, 
  Z, 
  b, 
  a, 
  t_miss = 0,
  s_miss = 0,
  sigma = NULL,
  include_residuals = TRUE
) {
  
  # Check missingness.
  if (t_miss + s_miss > 1) {
    stop("To ensure at least 1 outcome is always observed, the sum of
         target and surrogate missingness cannot exceed 1.0.")
  }
  
  # Default covariance matrix.
  if (is.null(sigma)) {
    sigma <- diag(2)
  }
  if (!is.matrix(sigma)) {
    stop("Covariance matrix required for sigma.")
  }
  
  # Ensure all structures are matrices.
  if (!is.matrix(X)) {X <- as.matrix(X)}
  if (!is.matrix(Z)) {Z <- as.matrix(Z)}
  b <- as.matrix(b, ncol = 1)
  a <- as.matrix(a, ncol = 1)
 
  # Observations.
  n <- nrow(X)
  
  # Linear predictors
  eta_t <- MMP(X, b)
  eta_s <- MMP(Z, a)
 
  # Residuals.
  if (include_residuals) {
    eig <- eigen(sigma, only.values = TRUE)$values
    if (min(eig) <= 0) {
      stop("Covariance matrix is not positive definite.")
    }
    
    et <- stats::rnorm(
      n = n,
      mean = 0,
      sd = sqrt(sigma[1, 1])
    )
    
    es <- stats::rnorm(
      n = n,
      mean = sigma[2, 1] / sigma[1, 1] * et,
      sd = sqrt(sigma[2, 2] - sigma[2, 1] / sigma[1, 1] * sigma[1, 2])
    )
    
    ts_resid <- cbind(et, es)
  } else {
    ts_resid <- array(0, dim = c(n, 2))
  }
  
  # Outcomes.
  y <- cbind(eta_t, eta_s) + ts_resid

  # Target missingness.
  nt_miss <- floor(t_miss * n)
  if (nt_miss > 0) {
    draw <- sort(sample(x = n, size = nt_miss, replace = FALSE))
    y[draw, 1] <- NA
  }
  
  # Surrogate missingness.
  ns_miss <- floor(s_miss * n)
  if (ns_miss > 0) {
    
    # Remove subjects with missing target outcome as candidates.
    candidates <- seq_len(n)[!is.na(y[, 1])]
    draw <- sort(sample(x = candidates, size = ns_miss, replace = FALSE))
    y[draw, 2] <- NA
  }

  # Output.
  dimnames(y) <- list(seq_len(n), c("Target", "Surrogate"))
  return(y)
}
