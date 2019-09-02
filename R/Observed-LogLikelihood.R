# Purpose: Function to calculate the observed data log likelihood
# Updated: 19/07/30

#' Observed Data Log Likelihood
#'
#' @param P List of partitioned data. See \code{\link{partSubj}}.
#' @param b Target regression coefficient.
#' @param a Surrogate regression coefficient
#' @param S Target-surrogate covariance
#'
#' @return Observed data log likelihood.
#' @export

obsLogLik = function(P,b,a,S){
  # Structure
  b = matrix(b,ncol=1);
  a = matrix(a,ncol=1);
  # Lambda
  L = matInv(S);
  n0 = P$Dims$n0;
  n1 = P$Dims$n1;
  n2 = P$Dims$n2;

  # Contribution of complete cases
  t0 = P$Complete$t0;
  s0 = P$Complete$s0;
  X0 = P$Complete$X0;
  Z0 = P$Complete$Z0;

  E0 = cbind(t0-MMP(X0,b),s0-MMP(Z0,a));
  l0 = n0*log(matDet(S))+tr(MMP(L,matIP(E0,E0)));

  # Contribution of subjects with target missingness
  s1 = P$TarMiss$s1;
  Z1 = P$TarMiss$Z1;

  E1 = matrix(c(s1-MMP(Z1,a)),ncol=1);
  l1 = n1*log(S[2,2])+as.numeric(matIP(E1,E1))/S[2,2];

  # Contribution of subjects with surrogate missingness
  t2 = P$SurMiss$t2;
  X2 = P$SurMiss$X2;

  E2 = matrix(c(t2-MMP(X2,b)),ncol=1);
  l2 = n2*log(S[1,1])+as.numeric(matIP(E2,E2))/S[1,1];

  # Final log likelihood
  l = -0.5*(l0+l1+l2);

  # Return
  return(l);
}
