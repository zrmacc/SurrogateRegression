# Purpose: Tabulate regression and covariance estimates and standard errors.
# Updated: 2022-08-05


#' Tabulate Regression Coefficients
#'
#' @param point Point estimates.
#' @param info Information matrix.
#' @param sig Significance level.
#' @return Data.table containing the point estimate, standard error, confidence
#'   interval, and Wald p-value.
RegTab <- function(point, info, sig = 0.05) {
  
  # Standard errors.
  se <- sqrt(diag(matInv(info)))
  
  # Critical value.
  z <- stats::qnorm(p = 1 - (sig / 2))
  
  # CIs.
  lower <- point - z * se
  upper <- point + z * se
  
  # P-values.
  p <- 2 * stats::pnorm(q = abs(point / se), lower.tail = FALSE)
  
  # Output.
  out <- data.frame(
    "Coefficient" = names(point), 
    "Point" = point, 
    "SE" = se, 
    "L" = lower, 
    "U" = upper, 
    "p" = p
  )
  rownames(out) <- seq_len(nrow(out))
  return(out)
}


# -----------------------------------------------------------------------------


#' Tabulate Covariance Parameters
#'
#' @param point Point estimates.
#' @param info Information matrix.
#' @param sig Significance level.
#' @return Data.table containing the point estimate, standard error, and
#'   confidence interval.
CovTab <- function(point, info, sig = 0.05) {
  
  # Critical value.
  z <- stats::qnorm(p = 1 - (sig / 2))
  
  # Inverse information.
  infoi <- matInv(info)
  
  # Standard errors
  se <- sqrt(diag(infoi))
  
  # Form CI for target-surrogate covariance on linear scale,
  # since this parameter is two-sided.
  lower2 <- point[2] - se[2]
  upper2 <- point[2] + se[2]

  # Log-scale standard errors for target-target and surrogate-surrogate
  # covariance, since these parameters are strictly positive.
  J <- diag(point[c(1, 3)])
  log.point <- log(point[c(1, 3)])
  log.info <- matQF(X = J, A = info[c(1, 3), c(1, 3)])
  log.infoi <- matInv(log.info)
  log.se <- sqrt(diag(log.infoi))

  # CI for variances on log scale.
  lower13 <- exp(log.point - z * log.se)
  upper13 <- exp(log.point + z * log.se)
  lower <- c(lower13[1], lower2, lower13[2])
  upper <- c(upper13[1], upper2, upper13[2])
  
  # Output
  out <- data.frame(
    "Covariance" = names(point), 
    "Point" = point, 
    "SE" = se, 
    "L" = lower, 
    "U" = upper
  )
  rownames(out) <- seq_len(nrow(out))
  return(out)
}


# -----------------------------------------------------------------------------


#' Format Output
#' 
#' @param data_part List of partitioned data. See \code{\link{PartitionData}}.
#' @param method Estimation method.
#' @param b Final target regression parameter.
#' @param a Final surrogate regression parameter.
#' @param sigma Final target-surrogate covariance matrix.
#' @param sig Significance level.
#' @return Object of class 'bnr'.
FormatOutput <- function(
  data_part,
  method,
  b,
  a,
  sigma,
  sig
) { 
  
  # Format sigma.
  colnames(sigma) <- rownames(sigma) <- c("Target", "Surrogate")
  
  # Covariance information.
  cov_info <- CovInfo(data_part = data_part, sigma = sigma)
  colnames(cov_info) <- 
    rownames(cov_info) <- 
      c("Target-Target", "Target-Surrogate", "Surrogate-Surrogate")
  
  # Covariance parameter table.
  cov_est <- c(
    "Target" = sigma[1, 1], 
    "Target-Surrogate" = sigma[1, 2], 
    "Surrogate" = sigma[2, 2]
  )
  cov_tab <- CovTab(point = cov_est, info = cov_info, sig = sig)
  
  # -----------------------------------
  
  # Regression information.
  reg_info <- RegInfo(data_part = data_part, sigma = sigma, as_matrix = TRUE)
  
  # Formatting.
  reg_est <- c(b, a)
  names(reg_est) <- c(data_part$Dims$x_names, data_part$Dims$z_names)
  
  # Regression coefficient table.
  reg_tab <- RegTab(reg_est, reg_info, sig = sig)
  reg_tab$Outcome <- c(rep("Target", data_part$Dims$q), rep("Surrogate", data_part$Dims$r))
  reg_tab <- reg_tab[, c(7, 1:6)]
  colnames(reg_info) <- rownames(reg_info) <- reg_tab$Coefficient
  
  # -----------------------------------
  
  # Residuals
  resid_mat <- cbind(
    data_part$Orig$t - MMP(data_part$Orig$X, matrix(b, ncol = 1)), 
    data_part$Orig$s - MMP(data_part$Orig$Z, matrix(a, ncol = 1))
  )
  colnames(resid_mat) <- c("Target", "Surrogate")
  rownames(resid_mat) <- seq_len(nrow(resid_mat))
  
  # Output
  out <- methods::new(
    Class = "bnr",
    Covariance = sigma,
    Covariance.info = cov_info,
    Covariance.tab = cov_tab,
    Method = method,
    Regression.info = reg_info,
    Regression.tab = reg_tab,
    Residuals = resid_mat
  )
  
  return(out)
}
