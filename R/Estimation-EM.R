# Purpose: Fitting procedure for bivariate normal regression via EM
# Updated: 2022-08-05

# -----------------------------------------------------------------------------
# Initialization.
# -----------------------------------------------------------------------------

#' Parameter Initialization
#'
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param b0 Initial target regression coefficient.
#' @param a0 Initial surrogate regression coefficient.
#' @param sigma0 Initial covariance matrix.
#' @return List containing initial values of beta, alpha, sigma.
ParamInit <- function(data_part, b0, a0, sigma0) {
 
  # Dimensions.
  q <- data_part$Dims$q
  r <- data_part$Dims$r
  n0 <- data_part$Dims$n0
  n1 <- data_part$Dims$n1
  n2 <- data_part$Dims$n2
  
  # -----------------------------------
  
  # Initialize beta.
  out <- list()
  if (is.null(b0)) {
    A <- array(0, dim = c(q, q))
    b <- array(0, dim = c(q, 1))
    
    # Contribution of complete cases.
    if (n0 > 0) {
      X0 <- data_part$Complete$X0
      t0 <- data_part$Complete$t0
      A <- A + matIP(X0, X0)
      b <- b + matIP(X0, t0)
    }
    
    # Contribution of surrogate missingness.
    if (n2 > 0) {
      X2 <- data_part$SMiss$X2
      t2 <- data_part$SMiss$t2
      A <- A + matIP(X2, X2)
      b <- b + matIP(X2, t2)
    }
    out$b <- solve(A, b)
    rm(A, b)
  } else {
    out$b <- b0
    rm(b0)
  }
  
  # -----------------------------------

  # Initialize alpha.
  if (is.null(a0)) {
    A <- array(0, dim = c(r, r))
    b <- array(0, dim = c(r, 1))
    
    # Contribution of complete cases.
    if (n0 > 0) {
      Z0 <- data_part$Complete$Z0
      s0 <- data_part$Complete$s0
      A <- A + matIP(Z0, Z0)
      b <- b + matIP(Z0, s0)
    }
    
    # Contribution of target missingness.
    if (n1 > 0) {
      Z1 <- data_part$TMiss$Z1
      s1 <- data_part$TMiss$s1
      A <- A + matIP(Z1, Z1)
      b <- b + matIP(Z1, s1)
    }

    out$a <- solve(A, b)
    rm(A, b)

  } else {

    out$a <- a0
    rm(a0)

  }
  
  # -----------------------------------

  # Initialize sigma.
  if (is.null(sigma0)) {
    
    # Only complete cases can contribute to initiation of sigma.
    if (n0 > 0) {
      E0 <- cbind(t0 - MMP(X0, out$b), s0 - MMP(Z0, out$a))
      out$sigma <- matIP(E0, E0) / n0
      rm(E0)
    } else {
  
      # If no complete cases are present, initialize to a diagonal matrix.
      out$sigma <- diag(
        stats::var(t0 - MMP(X0, out$b)),
        stats::var(s0 - MMP(Z0, out$a))
      )

    }
  } else {

    out$sigma <- sigma0
    rm(sigma0)

  }
  
  # Output.
  return(out)
}


# -----------------------------------------------------------------------------

#' EM Update
#' 
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param b0 Initial target regression coefficient.
#' @param a0 Initial surrogate regression coefficient.
#' @param sigma0 Initial covariance matrix.
#' 
#' @return List containing updated values for beta 'b', alpha 'a', 'sigma', the
#'   log likelihood 'loglik', and the change in log likelihood 'delta'.
UpdateEM <- function(
  data_part,
  b0,
  a0,
  sigma0
) {
  
  # Initial log likelihood.
  loglik0 <- ObsLogLik(
    data_part = data_part,
    b = b0, 
    a = a0, 
    sigma = sigma0
  )
  
  # Update regression parameters.
  out <- RegUpdate(data_part = data_part, sigma = sigma0)
  
  # Update covariance matrix.
  out$sigma <- CovUpdate(
    data_part = data_part,
    b0 = b0, 
    a0 = a0, 
    b1 = out$b, 
    a1 = out$a, 
    sigma0 = sigma0
  )
  
  # Final log likelihood
  out$loglik <- ObsLogLik(
    data_part = data_part,
    b = out$b, 
    a = out$a, 
    sigma = out$sigma
  )
  
  # Output
  out$delta <- (out$loglik - loglik0)
  return(out)
}


# -----------------------------------------------------------------------------

#' Update Iteration
#' 
#' @param theta0 List containing the initial parameter values.
#' @param update Function to iterate. Should accept and return a list similar
#'   parameter values.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report fitting progress?
IterUpdate <- function(
  theta0,
  update,
  maxit,
  eps,
  report
) {
  
  # Maximzation
  for (i in 1:maxit) {
    
    # Update.
    theta1 <- update(theta0)
    
    # Accept update if log likelihood has increased.
    if (theta1$delta > 0) {
      theta0 <- theta1
      if (report) {
        cat("Objective increment: ", signif(theta1$delta, digits = 3), "\n")
      }
    }
    
    # Terminate if increment is below tolerance
    if (theta1$delta < eps) {
      break
    }
  }
  
  # Fitting report
  if (report) {
    if (i < maxit) {
      cat(paste0(i - 1, " update(s) performed before tolerance limit.\n\n"))
    } else {
      cat(paste0(i, " update(s) performed without reaching tolerance limit.\n\n"))
    }
  }
  
  # Output.
  return(theta0)
}


# -----------------------------------------------------------------------------


#' Fit Bivariate Normal Regression Model via Expectation Maximization.
#'
#' Estimation procedure for bivariate normal regression models in which
#' the target and surrogate outcomes are both subject to missingness.
#'
#' The target and surrogate model matrices are expected in numeric format.
#' Include an intercept if required. Expand factors and interactions in advance.
#' Initial values may be specified for any of the target coefficient
#' \code{b0}, the surrogate coefficient \code{a0}, or the target-surrogate
#' covariance matrix \code{sigma0}.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @param b0 Initial target regression coefficient.
#' @param a0 Initial surrogate regression coefficient.
#' @param sigma0 Initial covariance matrix.
#' @param sig Type I error level.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report fitting progress?
#' @return An object of class 'bnr' with slots containing the estimated
#'   regression coefficients, the target-surrogate covariance matrix, the
#'   information matrices for the regression and covariance parameters, and the
#'   residuals.
FitBNEM <- function(
  t, 
  s, 
  X, 
  Z, 
  sig = 0.05, 
  b0 = NULL, 
  a0 = NULL, 
  sigma0 = NULL,
  maxit = 100, 
  eps = 1e-6, 
  report = TRUE
) {
  
  # Partition data.
  data_part <- PartitionData(t, s, X, Z)

  # Initialize parameters.
  theta0 <- ParamInit(data_part, b0, a0, sigma0)
  
  # Update wrapper.
  Update <- function(theta) {
    out <- UpdateEM(
      data_part = data_part,
      b0 = theta$b,
      a0 = theta$a,
      sigma0 = theta$sigma
    )
    return(out)
  }

  # Maximize.
  theta1 <- IterUpdate(
    theta0 = theta0,
    update = Update,
    maxit = maxit,
    eps = eps,
    report = report
  )
  
  # -----------------------------------

  # Output
  out <- FormatOutput(
    data_part = data_part,
    method = "EM",
    b = theta1$b,
    a = theta1$a,
    sigma = theta1$sigma,
    sig = sig
  )
  return(out)
}
