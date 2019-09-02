# Purpose: Wald test for bivariate normal regression via EM.
# Updated: 19/08/02

#' Wald Test via Expectation Maximization.
#'
#' Performs a Wald test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @param L Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param init Optional list of initial parameters for fitting the null model,
#'   with one or more of the components: a0, b0, S0.
#' @param maxit Maximum number of parameter updates.
#' @param eps Minimum acceptable improvement in log likelihood.
#' @param report Report model fitting progress? Default is FALSE.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#'
#' @return A numeric vector containing the Wald statistic, the degrees of
#'   freedom, and a p-value.

Wald.bnem = function(t,s,X,Z,L,init=NULL,maxit=100,eps=1e-8,report=F){
  # Input check
  if((!is.null(init))&&(!is.list(init))){stop("If initial parameter are provided, init should take the form
                                              of a list with one or more of the elements a0, b0, S0")};
  # Test specification
  p = ncol(X);
  df = sum(L);
  if(length(L)!=p){stop("L should have one entry per column of X.")};
  if(df==0){stop("At least 1 entry of L should be TRUE.")};
  if(df==p){stop("At least 1 entry of L should be FALSE.")};

  ## Model Fitting
  M0 = fit.bnem(t=t,s=s,X=X,Z=Z,b0=init$b0,a0=init$a0,S0=init$S0,maxit=maxit,eps=eps,report=report);

  # Extract information
  I = vcov(M0,type="Regression",inv=F);
  # Surrogate covariates
  q = ncol(Z);

  # Information keys
  key0 = c(L,rep(F,q));
  key1 = c(!L,rep(T,q));

  # Partition information
  Ibb = I[key0,key0,drop=F];
  Iaa = I[key1,key1,drop=F];
  Iba = I[key0,key1,drop=F];
  # Efficient information
  V = SchurC(Ibb=Ibb,Iaa=Iaa,Iba=Iba);

  ## Test
  # Coefficients of interest
  U = matrix(coef(M0,type="Target")[L,"Point"],ncol=1);
  # Statistic
  Tw = as.numeric(matQF(X=U,A=V));
  # P value
  p = pchisq(q=Tw,df=df,lower.tail=F);
  # Output
  Out = c("Wald"=Tw,"df"=df,"p"=p);
  return(Out);
}
