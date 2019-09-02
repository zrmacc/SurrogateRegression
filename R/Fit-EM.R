# Purpose: Fitting procedure for bivariate normal regression via EM
# Updated: 19/07/30

#' Parameter Initialization
#'
#' @param P List of partitioned data. See \code{\link{partSubj}}.
#' @param b0 Initial target regression coefficient.
#' @param a0 Initial surrogate regression coefficient.
#' @param S0 Initial covariance matrix.
#' @return List containing initial values of beta, alpha, sigma.

paramInit = function(P,b0,a0,S0){
  # Output structure
  Out = list();
  q = P$Dims$q;
  r = P$Dims$r;
  n0 = P$Dims$n0;
  n1 = P$Dims$n1;
  n2 = P$Dims$n2;

  # Beta
  if(is.null(b0)){
    A = array(0,dim=c(q,q));
    b = array(0,dim=c(q,1));
    if(n0>0){
      X0 = P$Complete$X0;
      t0 = P$Complete$t0;
      A = A+matIP(X0,X0);
      b = b+matIP(X0,t0);
    }
    if(n2>0){
      X2 = P$SurMiss$X2;
      t2 = P$SurMiss$t2;
      A = A+matIP(X2,X2);
      b = b+matIP(X2,t2);
    }
    Out$b = solve(A,b);
    rm(A,b);
  } else {
    Out$b = b0;
    rm(b0);
  };

  # Alpha
  if(is.null(a0)){
    A = array(0,dim=c(r,r));
    b = array(0,dim=c(r,1));
    if(n0>0){
      Z0 = P$Complete$Z0;
      s0 = P$Complete$s0;
      A = A+matIP(Z0,Z0);
      b = b+matIP(Z0,s0);
    }
    if(n1>0){
      Z1 = P$TarMiss$Z1;
      s1 = P$TarMiss$s1;
      A = A+matIP(Z1,Z1);
      b = b+matIP(Z1,s1);
    }
    Out$a = solve(A,b);
    rm(A,b);
  } else {
    Out$a = a0;
    rm(a0);
  }

  # Sigma
  if(is.null(S0)){
    if(n0>0){
      E0 = cbind(t0-MMP(X0,Out$b),s0-MMP(Z0,Out$a));
      Out$S = matIP(E0,E0)/n0;
      rm(E0);
    } else {
      warning("If no observations are complete, initializing the covariance is preferred.\n");
      Out$S = diag(2);
    }
  } else {
    Out$S = S0;
    rm(S0);
  }
  Out$L = matInv(Out$S);

  # Return
  return(Out)
}

#' Fit Bivariate Normal Regression Model via Expectation Maximization.
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
#' @return An object of class 'bnr' with slots containing the estimated
#'   regression coefficients, the target-surrogate covariance matrix, the
#'   information matrices for the regression and covariance parameters, and the
#'   residuals.
#' @examples
#' \dontrun{
#' # See `? rBNR` for data generation
#' fit.bnem(t=t,s=s,X=X,Z=Z,maxit=100,eps=1e-8,report=T);
#' }

fit.bnem = function(t,s,X,Z,sig=0.05,b0=NULL,a0=NULL,S0=NULL,maxit=100,eps=1e-6,report=T){
  # Partition subjects
  P = partSubj(t,s,X,Z);

  ## Initialize parameters
  theta0 = paramInit(P,b0,a0,S0);
  # Initial log likelihood
  theta0$ll = obsLogLik(P=P,b=theta0$b,a=theta0$a,S=theta0$S);

  ## Update Function
  Update = function(theta){
    # Update regression parameters
    Out = regUpdate(P=P,S=theta$S);

    # Update covariance matrix
    Out$S = covUpdate(P=P,b0=theta$b,a0=theta$a,b1=Out$b,a1=Out$a,S0=theta$S);

    # Final log likelihood
    Out$ll = obsLogLik(P=P,b=Out$b,a=Out$a,S=Out$S);

    ## Output
    Out$d = (Out$ll-theta$ll);
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

  ## Covariance information
  CovInfo = covInfo(P=P,S=S1);

  # Covariance parameter table
  Point = c(S1[1,1],S1[1,2],S1[2,2]);
  names(Point) = c("Target","Target-Surrogate","Surrogate");
  CovTab = covTab(point=Point,info=CovInfo,sig=sig);

  ## Regression information
  RegInfo = regInfo(P=P,S=S1);
  RegInfo = rbind(cbind(RegInfo$Ibb,RegInfo$Iba),cbind(t(RegInfo$Iba),RegInfo$Iaa));

  # Formatting
  Point = c(theta0$b,theta0$a);
  if(is.null(colnames(X))){colnames(X) = paste0("x",seq(1:ncol(X)))};
  if(is.null(colnames(Z))){colnames(Z) = paste0("z",seq(1:ncol(Z)))};
  names(Point) = c(colnames(X),colnames(Z));

  # Regression coefficient table
  RegTab = regTab(Point,RegInfo,sig=sig);
  RegTab$Outcome = c(rep("Target",P$Dims$q),rep("Surrogate",P$Dims$r));
  RegTab = RegTab[,c(7,1:6)];

  ## Residuals
  E = cbind(t-MMP(X,matrix(theta0$b,ncol=1)),s-MMP(Z,matrix(theta0$a,ncol=1)));
  colnames(E) = c("Target","Surrogate");

  ## Output
  Out = new(Class="bnr",
            Covariance=S1,
            Covariance.info=CovInfo,
            Covariance.tab=CovTab,
            Regression.info=RegInfo,
            Regression.tab=RegTab,
            Residuals=E);
  return(Out);
}
