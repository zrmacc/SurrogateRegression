# Purpose: Score test for bivariate normal regression via EM.
# Updated: 2020-11-28.

#' Score Test via Expectation Maximization.
#'
#' Performs a Score test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @param is_zero Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param init Optional list of initial parameters for fitting the null model.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#'
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value.
#'   
#' @examples 
#' \donttest{
#' # Generate data.
#' set.seed(100)
#' n <- 1e3
#' X <- cbind(1, rnorm(n))
#' Z <- cbind(1, rnorm(n))
#' data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.1, s_miss = 0.1)
#' 
#' # Test 1st coefficient.
#' score_test1 <- Score.BNEM(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(TRUE, FALSE)
#' )
#' 
#' # Test 2nd coefficient.
#' score_test2 <- Score.BNEM(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(FALSE, TRUE)
#' )
#' }

Score.BNEM <- function(
  t, 
  s, 
  X, 
  Z, 
  is_zero, 
  init = NULL, 
  maxit = 100, 
  eps = 1e-8, 
  report = FALSE
) {
  
  # Input check
  CheckInit(init = init)
  p <- ncol(X)
  CheckTestSpec(is_zero = is_zero, p = p)

  # Partition
  X.test <- X[, is_zero, drop = FALSE]
  X.null <- X[, !is_zero, drop = FALSE]
  
  # Fit null model
  fit0 <- Fit.BNEM(
    t = t, 
    s = s, 
    X = X.null, 
    Z = Z, 
    b0 = init$b0, 
    a0 = init$a0, 
    sigma0 = init$S0, 
    maxit = maxit, 
    eps = eps, 
    report = report
  )
  
  # -----------------------------------

  # Extract covariance.
  sigma <- vcov(fit0, type = "Outcome", inv = FALSE)
  sigma_inv <- matInv(sigma)
  
  # Extract residuals.
  t_resid <- matrix(resid(fit0, type = "Target"), ncol = 1)
  s_resid <- matrix(resid(fit0, type = "Surrogate"), ncol = 1)

  # Partition data.
  part_data_null <- PartitionData(t = t_resid, s = s_resid, X = X.null, Z = Z)
  part_data_test <- PartitionData(t = t_resid, s = s_resid, X = X.test, Z = Z)
  
  n0 <- part_data_null$Dims$n0
  n2 <- part_data_null$Dims$n2
  
  # -----------------------------------

  # Score statistic.
  df <- sum(is_zero)
  score <- array(0, dim = c(df, 1))
  
  # Contribution of complete cases.
  if (n0 > 0) {
    score <- score + sigma_inv[1, 1] * matIP(
      part_data_test$Complete$X0, 
      part_data_test$Complete$t0
    )
    score <- score + sigma_inv[1, 2] * matIP(
      part_data_test$Complete$X0, 
      part_data_test$Complete$s0
    )
  }
  
  # Contribution of surrogate missing.
  if (n2 > 0) {
    score <- score + matIP(
      part_data_test$SMiss$X2,
      part_data_test$SMiss$t2
    ) / sigma[1, 1]
  }
  
  # -----------------------------------

  # Target information.
  ibb <- array(0, dim = c(df, df))
  
  # Contribution of complete cases.
  if (n0 > 0) {
    ibb <- ibb + sigma_inv[1, 1] * matIP(
      part_data_test$Complete$X0,
      part_data_test$Complete$X0
    )
  }
  
  # Contribution of surrogate missing.
  if (n2 > 0) {
    ibb <- ibb + matIP(
      part_data_test$SMiss$X2, 
      part_data_test$SMiss$X2
    ) / sigma[1, 1]
  }

  # -----------------------------------
  
  # Surrogate information.
  iaa <- vcov(fit0, type = "Regression")
  
  # -----------------------------------

  # Cross information.
  iba <- array(0, dim = c(df, ncol(iaa)))
  
  # Contribution of complete cases.
  if (n0 > 0) {
    iba <- iba + cbind(
      sigma_inv[1, 1] * matIP(
        part_data_test$Complete$X0, 
        part_data_null$Complete$X0
      ),
      sigma_inv[1, 2] * matIP(
        part_data_test$Complete$X0,
        part_data_test$Complete$Z0
      )
    )
  }
  
  # Contribution of surrogate missing.
  if (n2 > 0) {
    iba[1:df, 1:(p - df)] <- iba[1:df, 1:(p - df)] + 
      matIP(
        part_data_test$SMiss$X2, 
        part_data_null$SMiss$X2
      ) / sigma[1, 1]
  }
  score_var <- SchurC(Ibb = ibb, Iaa = iaa, Iba = iba)
  
  # -----------------------------------

  # Test statistic.
  test_stat <- as.numeric(matQF(X = score, A = matInv(score_var)))
  
  # P-value.
  pval <- pchisq(q = test_stat, df = df, lower.tail = FALSE)
  
  # Output.
  out <- c(
    "Score" = test_stat, 
    "df" = df, 
    "p" = pval
  )
  return(out)
}
