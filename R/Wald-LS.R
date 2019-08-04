# Purpose: Wald test for bivariate normal regression via least squares.
# Updated: 19/08/02

#' Wald Test via Least Squares.
#'
#' Performs a Wald test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Model matrix.
#' @param L Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#'
#' @return A numeric vector containing the Wald statistic, the degrees of
#'   freedom, and a p-value.

Wald.bnls = function(t,s,X,L){
  # Test specification
  p = ncol(X);
  df = sum(L);
  if(length(L)!=p){stop("L should have one entry per column of X.")};
  if(df==0){stop("At least 1 entry of L should be TRUE.")};
  if(df==p){stop("At least 1 entry of L should be FALSE.")};

  ## Model Fitting
  M0 = fit.bnls(t=t,s=s,X=X);

  # Extract information
  I = vcov(M0,type="Regression",inv=F);

  # Information keys
  key0 = c(L,rep(F,p));
  key1 = c(!L,rep(T,p));

  # Partition information
  Ibb = I[key0,key0,drop=F];
  Iaa = I[key1,key1,drop=F];
  Iba = I[key0,key1,drop=F];
  # Efficient information
  V = SchurC(Ibb=Ibb,Iaa=Iaa,Iba=Iba);

  ## Test
  # Coefficients of interest
  U = coef(M0,type="Target")[L,"Point"];
  # Statistic
  Tw = as.numeric(matQF(X=U,A=V));
  # P value
  p = pchisq(q=Tw,df=df,lower.tail=F);
  # Output
  Out = c("Wald"=Tw,"df"=df,"p"=p);
  return(Out);
}
