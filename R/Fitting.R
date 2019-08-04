# Purpose: Master fitting function for bivariate normal regression
# Updated: 19/06/18

#' Fit Bivariate Normal Regression Model.
#'
#' Estimation procedure for bivariate normal regression models. The EM algorithm
#' is applied if \code{s} contains missing values, or if \code{X} differs from
#' \code{Z}. Otherwise, an accelerated least squares procedure is applied.
#'
#' The target and surrogate model matrices are expected in numeric format.
#' Include an intercept if required. Expand factors and interactions in advance.
#' See \code{\link{fit.bnem}} for additional options accepted when fitting via EM.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix. Defaults to \code{X}.
#' @param sig Significance level.
#' @param ... Additional arguments accepted if fitting via EM.
#'
#' @importFrom methods new
#' @importFrom stats coef pnorm qnorm resid
#' @export
#' @return An object of class 'mnr' with slots containing the estimated regression
#'  coefficients, the target-surrogate covariance matrix, the information matrices
#'  for regression parameters, and the residuals.
#' @examples
#' \dontrun{
#' # See `? rBNR` for data generation
#' fit.bnr(t=t,s=s,X=X,Z=Z,maxit=100,eps=1e-8,report=T);
#' }

fit.bnr = function(t,s,X,Z=NULL,sig=0.05,...){
  # Input check
  if(!is.vector(t)){stop("A numeric vector is expected for t.")};
  if(!is.vector(s)){stop("A numeric vector is expected for s.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  if((!is.null(Z))&(!is.matrix(Z))){stop("A numeric matrix is expected for Z.")};

  # Determine if s contains missing values, or if Z differs from X.
  EM = (sum(is.na(s))>0)|((!is.null(Z))&(!identical(X,Z)));
  # If missingness occurs in s, apply EM algorithm.
  if(EM){
    if(is.null(Z)){Z=X};
    Out = fit.bnem(t=t,s=s,X=X,Z=Z,sig=sig,...);
  } else {
    # Otherwise, apply the least squares procedure.
    Out = fit.bnls(t=t,s=s,X=X);
  }
  # Output
  return(Out);
};
