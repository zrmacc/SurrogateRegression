# Purpose: ECM update for covariance parameters
# Updated: 19/08/01

#' Covariate Update
#'
#' @param P List of partitioned data. See \code{\link{partSubj}}.
#' @param b0 Previous target regression coefficient.
#' @param a0 Previous surrogate regression coefficient.
#' @param b1 Current target regression coefficient.
#' @param a1 Current surrogate regression coefficient.
#' @param S0 Initial target-surrogate covariance matrix.
#' @export
#' @return ECM update of the target-surrogate covariance matrix.

covUpdate = function(P,b0,a0,b1,a1,S0){
  ## Dimensions
  n0 = P$Dims$n0;
  n1 = P$Dims$n1;
  n2 = P$Dims$n2;

  # Lambda
  L0 = matInv(S0);

  ## Residual vectors
  et = es = c();
  # Complete cases
  if(n0>0){
    # Target
    et0 = P$Complete$t0-MMP(P$Complete$X0,b1);
    et = c(et,et0);
    # Surrogate
    es0 = P$Complete$s0-MMP(P$Complete$Z0,a1);
    es = c(es,es0);
  }
  # Target missing
  if(n1>0){
    # Surrogate
    es1 = P$TarMiss$s1-MMP(P$TarMiss$Z1,a1);
    es = c(es,es1);
    # Target
    w1 = (S0[1,2]/S0[2,2]);
    et1 = MMP(P$TarMiss$X1,b0-b1)+w1*(P$TarMiss$s1-MMP(P$TarMiss$Z1,a0));
    et = c(et,et1);
  }
  if(n2>0){
    # Target
    et2 = P$SurMiss$t2-MMP(P$SurMiss$X2,b1);
    et = c(et,et2);
    # Surrogate
    w2 = (S0[2,1]/S0[1,1]);
    es2 = MMP(P$SurMiss$Z2,a0-a1)+w2*(P$SurMiss$t2-MMP(P$SurMiss$X2,b0));
    es = c(es,es2);
  }

  # Residual matrix
  E = cbind(et,es);
  # Expected outer product
  V = matIP(E,E)+n1*diag(c(1/L0[1,1],0))+n2*diag(c(0,1/L0[2,2]));

  ## Update sigma
  S = (V/P$Dims$n);
  return(S);
}

