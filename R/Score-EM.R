# Purpose: Score test for bivariate normal regression via EM.
# Updated: 19/08/02

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
  X.test = X[,L,drop=F];
  X.null = X[,!L,drop=F];
  # Fit null model
  M0 = fit.bnem(t=t,s=s,X=X.null,Z=Z,b0=init$b0,a0=init$a0,S0=init$S0,maxit=maxit,eps=eps,report=report);

  # Extract covariance
  S = vcov(M0,type="Outcome",inv=F);
  L = matInv(S);
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = resid(M0,type="Surrogate");

  ## Partition
  P = list();
  P$Inds$TarMiss = is.na(eT);
  P$Inds$SurMiss = is.na(eS);
  n0 = sum(P$Inds$TarMiss);
  n2 = sum(P$Inds$SurMiss);

  # Keys
  P$Complete$key = (!P$Inds$TarMiss)&(!P$Inds$SurMiss);
  P$SurMiss$key = (!P$Inds$TarMiss)&(P$Inds$SurMiss);

  # Complete cases
  P$Complete$eT = eT[P$Complete$key];
  P$Complete$eS = eS[P$Complete$key];
  P$Complete$X.test = X.test[P$Complete$key,,drop=F];
  P$Complete$X.null = X.null[P$Complete$key,,drop=F];
  P$Complete$Z = Z[P$Complete$key,,drop=F];

  # Surrogate missing
  P$SurMiss$eT = eT[P$SurMiss$key];
  P$SurMiss$X.test = X.test[P$SurMiss$key,,drop=F];
  P$SurMiss$X.null = X.null[P$SurMiss$key,,drop=F];

  ## Score
  U = array(0,dim=c(df,1));
  if(n0>0){
    U = U + L[1,1]*matIP(P$Complete$X.test,P$Complete$eT);
    U = U + L[1,2]*matIP(P$Complete$X.test,P$Complete$eS);
  }
  if(n2>0){
    U = U + matIP(P$SurMiss$X.test,P$SurMiss$eT)/S[1,1];
  }

  ## Information
  # Target information
  Ibb = array(0,dim=c(df,df));
  if(n0>0){
    Ibb = Ibb + L[1,1]*matIP(P$Complete$X.test,P$Complete$X.test);
  }
  if(n2>0){
    Ibb = Ibb + matIP(P$SurMiss$X.test,P$SurMiss$X.test)/S[1,1];
  }

  # Nuisance information
  Iaa = vcov(M0,type="Regression");

  # Cross information
  Iba = array(0,dim=c(df,ncol(Iaa)));
  if(n0>0){
    Iba = Iba + cbind(L[1,1]*matIP(P$Complete$X.test,P$Complete$X.null),
                      L[1,2]*matIP(P$Complete$X.test,P$Complete$Z));
  }
  if(n2>0){
    Iba[1:df,1:(p-df)] = Iba[1:df,1:(p-df)] + matIP(P$SurMiss$X.test,P$SurMiss$X.null)/S[1,1];
  }
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
