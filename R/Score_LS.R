# Purpose: Score test for bivariate normal regression via least squares.
# Updated: 19/06/18

#' Score Test via Least Squares.
#'
#' Performs a Score test of the null hypothesis that a subset of the regression
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
#' @return A numeric vector containing the score statistic, the degrees of
#'   freedom, and a p-value.

Score.bnls = function(t,s,X,L){
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
  M0 = fit.bnls(t=t,s=s,X=Xb);
  # Extract covariance
  S = vcov(M0,type="Outcome",inv=F);
  L = matInv(S);
  # Extract residuals
  eT = resid(M0,type="Target");
  eS = resid(M0,type="Surrogate");

  ## Keys
  mt = is.na(eT);
  # Complete cases
  key0 = (!mt);

  # Complete cases
  eT0 = eT[key0];
  eS0 = eS[key0];
  Xa0 = Xa[key0,,drop=F];
  Xb0 = Xb[key0,,drop=F];

  ## Score
  U = array(0,dim=c(df,1));
  U = U+L[1,1]*matIP(Xa0,eT0)+L[1,2]*matIP(Xa0,eS0);

  ## Information
  # Target information
  Ibb = L[1,1]*matIP(Xa0,Xa0);
  # Nuisance information
  Iaa = vcov(M0,type="Information");
  # Cross information
  Iba = cbind(L[1,1]*matIP(Xa0,Xb0),L[1,2]*matIP(Xa0,Xb0));
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
