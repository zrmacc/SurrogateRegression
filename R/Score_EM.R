# Purpose: Score test for bivariate normal regression via EM.
# Updated: 18/09/11

#' Score Test via Expectation Maximization.
#'
#' Performs a Score test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @param L Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param init Optional list of initial parameters for fitting the null model.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#'
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value.

Score.bnem = function(t,s,X,Z,L,init=NULL,maxit=100,eps=1e-8,report=F){
  # Input check
  if((!is.null(init))&&(!is.list(init))){stop("If initial parameter are provided, init should take the form
                                               of a list with one or more of the elements a0, b0, S0")};
  # Test specification
  p = ncol(X);
  df = sum(L);
  if(length(L)!=p){stop("L should have one entry per column of X.")};
  if(df==0){stop("At least 1 entry of L should be TRUE.")};
  if(df==p){stop("At least 1 entry of L should be FALSE.")};

  ## Partition
  Xa = X[,L,drop=F];
  Xb = X[,!L,drop=F];
  # Fit null model
  M0 = fit.bnem(t=t,s=s,X=Xb,Z=Z,b0=init$b0,a0=init$a0,S0=init$S0,maxit=maxit,eps=eps,report=report);
  # Extract covariance
  S = vcov(M0,type="Outcome",inv=F);
  L = matInv(S);
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = resid(M0,type="Surrogate");

  ## Keys
  mt = is.na(eT);
  ms = is.na(eS);
  # Complete cases
  key0 = (!mt)&(!ms);
  # Surrogate missing
  key2 = (!mt)&(ms);

  # Complete cases
  eT0 = eT[key0];
  eS0 = eS[key0];
  Xa0 = Xa[key0,,drop=F];
  Xb0 = Xb[key0,,drop=F];
  Z0 = Z[key0,,drop=F];

  # Surrogate missing
  eT2 = eT[key2];
  Xa2 = Xa[key2,,drop=F];
  Xb2 = Xb[key2,,drop=F];

  ## Score
  U = array(0,dim=c(df,1));
  U = U+L[1,1]*matIP(Xa0,eT0)+L[1,2]*matIP(Xa0,eS0)+matIP(Xa2,eT2)/S[1,1];

  ## Information
  # Target information
  Ibb = L[1,1]*matIP(Xa0,Xa0)+matIP(Xa2,Xa2)/S[1,1];
  # Nuisance information
  Iaa = vcov(M0,type="Information");
  # Cross information
  Iba = cbind(L[1,1]*matIP(Xa0,Xb0)+matIP(Xa2,Xb2)/S[1,1],L[1,2]*matIP(Xa0,Z0));
  V = SchurC(Ibb=Ibb,Iaa=Iaa,Iba=Iba);

  ## Test
  # Statistic
  Ts = as.numeric(matQF(X=U,A=matInv(V)));
  # P value
  p = pchisq(q=Ts,df=df,lower.tail=F);
  # Output
  Out = c("Score"=Ts,"df"=df,"p"=p);
  return(Out);
}
