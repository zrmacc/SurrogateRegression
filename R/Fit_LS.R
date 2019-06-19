# Purpose: Fitting procedure for bivariate normal regression via LS
# Updated: 19/01/28

#' Fit Bivariate Outcome Model via Least Squares
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
#' @return An object of class 'bnr' with slots containing the estimated regression
#'  coefficients, the target-surrogate covariance matrix, the information matrices
#'  for regression parameters, and the residuals.
#' @examples
#' \dontrun{
#' # See `? rBNR` for data generation
#' fit.bnls(t=t,s=s,X=X,Z=Z,maxit=100,eps=1e-8,report=T);
#' } 

fit.bnls = function(t,s,X,sig=0.05){
  # Originals
  X00 = X;
  Y00 = cbind(t,s);
  
  # Covariate dimension
  p = ncol(X);
  # ID
  id = seq(1:length(t));
  
  ## Partition data
  mt = is.na(t);
  ms = is.na(s);
  # Check surrogate missingness
  if(sum(ms)>0){stop("No missingness in the target outcome is expected for BNLS.")};
  
  # Complete observations
  key0 = !(mt);
  t0 = t[key0];
  s0 = s[key0];
  X0 = X[key0,,drop=F];
  n0 = length(t0);
  id0 = id[key0];
  
  # Missing target and not surrogate
  key1 = mt;
  s1 = s[key1];
  X1 = X[key1,,drop=F];
  n1 = length(s1);
  id1 = id[key1];
  
  # Final sample size
  n = n0+n1;
  # Restructure
  X = rbind(X0,X1);
  s = c(s0,s1);
  
  # Stage 1 regression
  m1 = fitOLS(y=s,X=X);
  # Parameters
  alpha = m1$Beta;
  SigmaSS = m1$V;
  rm(m1);
  
  # Stage 2 regression
  Z0 = cbind(s0,X0);
  m2 = fitOLS(y=t0,X=Z0);
  zeta = m2$Beta;
  LambdaTTi = m2$V;
  rm(m2);
  
  # Recover original parameters
  delta = zeta[1];
  gamma = zeta[2:(p+1)];
  beta  = delta*alpha+gamma;
  SigmaTS = delta*SigmaSS;
  SigmaTT = LambdaTTi+SigmaSS*(delta)^2;
  LambdaSSi = SigmaSS-(SigmaTS^2)/SigmaTT;
  
  # Final covariance
  Sigma = matrix(c(SigmaTT,SigmaTS,SigmaTS,SigmaSS),nrow=2);
  colnames(Sigma) = rownames(Sigma) = c("Target","Surrogate");
  
  ## Information
  J = list();
  
  # For alpha
  J$Iaa = matIP(X0,X0)/LambdaSSi+matIP(X1,X1)/SigmaSS;
  J$Ibb = matIP(X0,X0)/LambdaTTi;
  J$Iba = -delta*J$Ibb;
  # Note: -(delta/LambdaTTi) = LambdaTS;
  
  # Overall
  Info = rbind(cbind(J$Ibb,J$Iba),cbind(t(J$Iba),J$Iaa));
  Infoi = matInv(Info);
  
  # Regression coefficients
  Point = c(beta,alpha);
  # SE
  SE = sqrt(diag(Infoi));
  # CIs
  z = qnorm(p=1-(sig/2));
  L = Point-z*SE;
  U = Point+z*SE;
  P = 2*pnorm(q=abs(Point/SE),lower.tail=F);
  # Labeling
  Outcome = c(rep("Target",p),rep("Surrogate",p));
  if(is.null(colnames(X))){colnames(X) = paste0("x",seq(1:ncol(X)))};
  Regressor = rep(colnames(X),times=2);
  Coeff = data.frame("Outcome"=Outcome,"Regressor"=Regressor,"Point"=Point,"SE"=SE,"L"=L,"U"=U,"p"=P);
  colnames(Info) = rownames(Info) = Coeff$Regressor;
  
  ## Residuals
  E = cbind(Y00[,1]-MMP(X00,beta),Y00[,2]-MMP(X00,alpha));
  colnames(E) = c("Target","Surrogate");
  
  ## Output
  Out = new(Class="bnr", Coefficients=Coeff,Covariance=Sigma,Information=Info,Residuals=E);
  return(Out);
}