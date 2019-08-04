# Purpose: Calculate information for covariance parameters
# Updated: 19/08/01

#' Covariance Information
#'
#' @param P List of partitioned data. See \code{\link{partSubj}}.
#' @param S Target-surrogate covariance matrix.
#' @return Information matrix for the residual target variance (Itt),
#'   target-surrogate covariance, and surrogate variance.

covInfo = function(P,S){
  # Dimensions
  n0 = P$Dims$n0;
  n1 = P$Dims$n1;
  n2 = P$Dims$n2;

  # Lambda
  L = matInv(S);

  Out = array(0,dim=c(3,3));
  # Complete cases
  if(n0>0){
    J = array(0,dim=c(3,3));
    J[1,1] = L[1,1]^2;
    J[2,2] = 2*(L[1,2]^2+L[1,1]*L[2,2]);
    J[3,3] = L[2,2]^2;
    J[1,2] = J[2,1] = 2*L[1,1]*L[1,2];
    J[2,3] = J[3,2] = 2*L[1,2]*L[2,2];
    J[1,3] = J[3,1] = L[1,2]^2;
    Out = Out + 0.5*n0*J;
  }
  if(n1>0){
    Out[3,3] = Out[3,3] + 0.5*n1/(S[2,2]^2);
  }
  if(n2>0){
    Out[1,1] = Out[1,1] + 0.5*n2/(S[1,1]^2);
  }
  return(Out);
}

