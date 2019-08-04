# Purpose: Fitting procedure for bivariate normal regression via LS
# Updated: 19/01/28

#' Fit Bivariate Normal Regression Model via Least Squares
#'
#' Estimation procedure for bivariate normal regression models in which
#' only the target outcome is subject to missingness.
#'
#' The model matrix is expected in numeric format. Include an intercept if
#' required. Expand factors and interactions in advance.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Model matrix.
#' @param sig Significance level.
#'
#' @importFrom methods new
#' @importFrom stats coef pnorm qnorm resid
#' @return An object of class 'bnr' with slots containing the estimated
#'   regression coefficients, the target-surrogate covariance matrix, the
#'   information matrices for the regression and covariance parameters, and the
#'   residuals.
#' @examples
#' \dontrun{
#' # See `? rBNR` for data generation
#' fit.bnls(t=t,s=s,X=X);
#' }

fit.bnls = function(t,s,X,sig=0.05){
  # Partition subjects
  P = partSubj(t,s,X);
  q = P$Dims$q;

  # Stage 1 regression
  fit.1 = fitOLS(y=s,X=X);

  # Parameters
  alpha = fit.1$Beta;
  S22 = fit.1$V;
  rm(fit.1);

  # Stage 2 regression
  Z0 = cbind(P$Complete$s0,P$Complete$X0);
  fit.2 = fitOLS(y=P$Complete$t0,X=Z0);
  zeta = fit.2$Beta;
  L11i = fit.2$V;
  rm(fit.2);

  # Recover original parameters
  delta = zeta[1];
  gamma = zeta[2:(q+1)];
  beta  = delta*alpha+gamma;
  S12 = delta*S22;
  S11 = L11i+S22*(delta)^2;

  # Final covariance
  S = matrix(c(S11,S12,S12,S22),nrow=2);
  colnames(S) = rownames(S) = c("Target","Surrogate");

  ## Covariance information
  CovInfo = covInfo(P=P,S=S);

  # Covariance parameter table
  Point = c(S[1,1],S[1,2],S[2,2]);
  names(Point) = c("Target","Target-Surrogate","Surrogate");
  CovTab = covTab(point=Point,info=CovInfo,sig=sig);

  ## Regression information
  RegInfo = regInfo(P=P,S=S);
  RegInfo = rbind(cbind(RegInfo$Ibb,RegInfo$Iba),cbind(t(RegInfo$Iba),RegInfo$Iaa));

  # Formatting
  Point = c(beta,alpha);
  if(is.null(colnames(X))){colnames(X) = paste0("x",seq(1:ncol(X)))};
  names(Point) = rep(colnames(X),times=2);

  # Regression coefficient table
  RegTab = regTab(Point,RegInfo,sig=sig);
  RegTab$Outcome = c(rep("Target",P$Dims$q),rep("Surrogate",P$Dims$r));
  RegTab = RegTab[,c(7,1:6)];

  colnames(RegInfo) = rownames(RegInfo) = RegTab$Regressor;

  ## Residuals
  E = cbind(t-MMP(X,beta),s-MMP(X,alpha));
  colnames(E) = c("Target","Surrogate");

  ## Output
  Out = new(Class="bnr",
            Covariance=S,
            Covariance.info=CovInfo,
            Covariance.tab=CovTab,
            Regression.info=RegInfo,
            Regression.tab=RegTab,
            Residuals=E);
  return(Out);
}
