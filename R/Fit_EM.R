# Purpose: Fitting procedure for bivariate normal regression via EM
# Updated: 19/01/28

#' Fit Bivariate Outcome Model via Expectation Maximization. 
#' 
#' Estimation procedure for bivariate normal regression models in which
#' the target and surrogate outcomes are both subject to missingness. 
#' 
#' The target and surrogate model matrices are expected in numeric format. 
#' Include an intercept if required. Expand factors and interactions in advance.
#' Initial values may be specified for any of the target coefficient
#' \code{b0}, the surrogate coefficient \code{a0}, or the target-surrogate
#' covariance matrix \code{S0}.
#'
#' @param t Target outcome vector. 
#' @param s Surrogate outcome vector. 
#' @param X Target model matrix. 
#' @param Z Surrogate model matrix. 
#' @param b0 Initial target regression coefficient.
#' @param a0 Initial surrogate regression coefficient.
#' @param S0 Initial covariance matrix.
#' @param sig Significance level. 
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report fitting progress? 
#'
#' @importFrom methods new
#' @importFrom stats coef pnorm qnorm resid
#' @return An object of class 'bnr' with slots containing the estimated regression
#'  coefficients, the target-surrogate covariance matrix, the information matrices
#'  for regression parameters, and the residuals.
#' @examples
#' \dontrun{
#' # See `? rBNR` for data generation
#' fit.bnem(t=t,s=s,X=X,Z=Z,maxit=100,eps=1e-8,report=T);
#' } 

fit.bnem = function(t,s,X,Z,sig=0.05,b0=NULL,a0=NULL,S0=NULL,maxit=100,eps=1e-6,report=T){
  # Originals
  X00 = X;
  Z00 = Z;
  # Covariate dimensions
  p = ncol(X);
  q = ncol(Z);
  # Dimensional consistency
  if(!is.null(b0)){
    if(length(b0)!=p){stop("A numeric vector of length ncol(X) is expected for b0.")};
  }
  if(!is.null(a0)){
    if(length(a0)!=q){stop("A numeric vector of length ncol(Z) is expected for a0.")};
  }
  if(!is.null(S0)){
    if(!is.matrix(S0)){stop("A numeric 2x2 target-surrogate covariance matrix is expected for S0.")};
  }

  # ID
  id = seq(1:length(t));
  ## Partition data
  mt = is.na(t);
  ms = is.na(s);
  # Complete observations
  flag = !(mt|ms);
  t0 = t[flag];
  s0 = s[flag];
  X0 = X[flag,,drop=F];
  Z0 = Z[flag,,drop=F];
  n0 = length(t0);
  id0 = id[flag];
  # Missing target and not surrogate
  flag = (!ms)&mt;
  s1 = s[flag];
  X1 = X[flag,,drop=F];
  Z1 = Z[flag,,drop=F];
  n1 = length(s1);
  id1 = id[flag];
  # Missing surrogate and not target
  flag = ms&(!mt);
  t2 = t[flag];
  X2 = X[flag,,drop=F];
  Z2 = Z[flag,,drop=F];
  n2 = length(t2);
  id2 = id[flag];
  # Missing both outcomes
  flag = ms&mt;
  id3 = id[flag];
  # Final sample size
  n = n0+n1+n2;
  # Restructure
  X = rbind(X0,X1,X2);
  Z = rbind(Z0,Z1,Z2);
  # Reusable products
  XtXi = matInv(matIP(X,X));
  ZtZi = matInv(matIP(Z,Z));
  
  ## Initialize
  theta0 = list();
  # Beta
  if(is.null(b0)){
    A = array(0,dim=c(p,p));
    b = array(0,dim=c(p,1));
    if(n0>0){
      A = A+matIP(X0,X0);
      b = b+matIP(X0,t0);
    }
    if(n2>0){
      A = A+matIP(X2,X2);
      b = b+matIP(X2,t2);
    }
    theta0$b = solve(A,b);
    rm(A,b);
  } else {
    theta0$b = b0;
    rm(b0);
  };
  
  # Alpha
  if(is.null(a0)){
    A = array(0,dim=c(q,q));
    b = array(0,dim=c(q,1));
    if(n0>0){
      A = A+matIP(Z0,Z0);
      b = b+matIP(Z0,s0);
    }
    if(n1>0){
      A = A+matIP(Z1,Z1);
      b = b+matIP(Z1,s1);
    }
    theta0$a = solve(A,b);
    rm(A,b);
  } else {
    theta0$a = a0;
    rm(a0);
  }
  
  # Sigma
  if(is.null(S0)){
    if(n0>0){
      E0 = cbind(t0-MMP(X0,theta0$b),s0-MMP(Z0,theta0$a));
      theta0$S = matIP(E0,E0)/n0;
      rm(E0);
    } else {
      warning("If no observations are complete, initializing the covariance is preferred.\n");
      theta0$S = diag(2);
    }
  } else {
    theta0$S = S0;
    rm(S0);
  } 
  theta0$L = matInv(theta0$S);
  
  ## Update Function
  Update = function(theta){
    # Current covariance
    S0 = theta$S;
    L0 = theta$L;
    
    ## Working vectors
    t = s = c();
    # Complete cases
    if(n0>0){
      # Surrogate
      s = c(s,s0);
      # Target
      t = c(t,t0);
    }
    # Target missing
    if(n1>0){
      # Surrogate
      s = c(s,s1);
      # Target
      w1 = (S0[1,2]/S0[2,2]);
      t1 = MMP(X1,theta$b)+w1*(s1-MMP(Z1,theta$a));
      t = c(t,t1);
    }
    if(n2>0){
      # Surrogate
      w2 = (S0[2,1]/S0[1,1]);
      s2 = MMP(Z2,theta$a)+w2*(t2-MMP(X2,theta$b));
      s = c(s,s2);
      # Target
      t = c(t,t2);
    }

    # Surrogate residual
    es = s-MMP(Z,theta$a);
    # Residual matrix
    E0 = cbind(t-MMP(X,theta$b),es);
    # Expected outer product
    V0 = matIP(E0,E0)+n1*diag(c(1/L0[1,1],0))+n2*diag(c(0,1/L0[2,2]));
    # Baseline objective
    Q0 = -n*log(det(S0))-tr(MMP(L0,V0));
    
    ## Update beta
    w1 = L0[1,2]/L0[1,1];
    b1 = MMP(XtXi,matIP(X,t+w1*es));
    
    # Target residual
    et = t-MMP(X,b1);
    
    ## Update alpha
    w2 = L0[1,2]/L0[2,2];
    a1 = MMP(ZtZi,matIP(Z,w2*et+s));
    
    # Final residual matrix
    E1 = cbind(et,s-MMP(Z,a1));
    # Expected outer product
    V1 = matIP(E1,E1)+n1*diag(c(1/L0[1,1],0))+n2*diag(c(0,1/L0[2,2]));
    
    ## Update sigma
    S1 = (V1/n);
    L1 = matInv(S1);
    
    # Final objective
    Q1 = -n*log(det(S1))-tr(MMP(L1,V1));
    d = Q1-Q0;
    
    ## Output
    Out = list("b"=b1,"a"=a1,"S"=S1,"L"=L1,"d"=d);
    return(Out);
  }
  
  ## Maximzation
  for(i in 1:maxit){
    # Update
    theta1 = Update(theta0);
    # Accept if increment is positive
    if(theta1$d>0){
      theta0 = theta1;
      if(report){cat("Objective increment: ",signif(theta1$d,digits=3),"\n")}
    }
    # Terminate if increment is below tolerance
    if(theta1$d<eps){
      rm(theta1);
      break;
    };
  };
  
  ## Fitting report
  if(report){
    if(i<maxit){
      cat(paste0(i-1," update(s) performed before tolerance limit.\n\n"));
    } else {
      cat(paste0(i," update(s) performed without reaching tolerance limit.\n\n"));
    };
  };
  
  # Final covariance
  S1 = theta0$S;
  colnames(S1) = rownames(S1) = c("Target","Surrogate");
  # Final precision
  L1 = theta0$L;
  
  ## Information
  J = list();
  # For Beta
  Ibb = array(0,dim=c(p,p));
  if(n0>0){
    Ibb = Ibb+L1[1,1]*matIP(X0,X0);
  }
  if(n2>0){
    Ibb = Ibb+matIP(X2,X2)/S1[1,1];
  }
  J$Ibb = Ibb;
  # For alpha
  Iaa = array(0,dim=c(q,q));
  if(n0>0){
    Iaa = Iaa+L1[2,2]*matIP(Z0,Z0);
  }
  if(n1>0){
    Iaa = Iaa+matIP(Z1,Z1)/S1[2,2];
  }
  J$Iaa = Iaa;
  # Cross information
  Iba = array(0,dim=c(p,q));
  if(n0>0){
    Iba = Iba+L1[1,2]*matIP(X0,Z0);
  }
  J$Iba = Iba;
  # Overall
  Info = rbind(cbind(J$Ibb,J$Iba),cbind(t(J$Iba),J$Iaa));
  Infoi = matInv(Info);
  
  # Regression coefficients
  Point = c(theta0$b,theta0$a);
  # SE
  SE = sqrt(diag(Infoi));
  # CIs
  z = qnorm(p=1-(sig/2));
  L = Point-z*SE;
  U = Point+z*SE;
  P = 2*pnorm(q=abs(Point/SE),lower.tail=F);
  # Labeling
  Outcome = c(rep("Target",p),rep("Surrogate",q));
  if(is.null(colnames(X))){colnames(X) = paste0("x",seq(1:ncol(X)))};
  if(is.null(colnames(Z))){colnames(Z) = paste0("z",seq(1:ncol(Z)))}; 
  Regressor = c(colnames(X),colnames(Z));
  Coeff = data.frame("Outcome"=Outcome,"Regressor"=Regressor,"Point"=Point,"SE"=SE,"L"=L,"U"=U,"p"=P);
  colnames(Info) = rownames(Info) = Coeff$Regressor;
  
  ## Residuals
  E = cbind(t-MMP(X00,theta0$b),s-MMP(Z00,theta0$a));
  colnames(E) = c("Target","Surrogate");
  
  ## Output
  Out = new(Class="bnr", Coefficients=Coeff,Covariance=S1,Information=Info,Residuals=E);
  return(Out);
}