# Purpose: Calculate information for regression parameters
# Updated: 19/08/01

#' Regression Information
#'
#' @param P List of partitioned data. See \code{\link{partSubj}}.
#' @param S Target-surrogate covariance matrix.
#' @return List containing the information matrix for beta (Ibb), the
#'   information matrix for alpha (Iaa), and the cross information (Iba).

regInfo = function(P,S){
  # Dimensions
  p = P$Dims$p;
  q = P$Dims$q;
  r = P$Dims$r;
  n0 = P$Dims$n0;
  n1 = P$Dims$n1;
  n2 = P$Dims$n2;

  # Lambda
  L = matInv(S);

  Out = list();
  # Information for beta
  Ibb = array(0,dim=c(q,q));
  if(n0>0){
    Ibb = Ibb+L[1,1]*P$IPs$X0tX0;
  }
  if(n2>0){
    Ibb = Ibb+P$IPs$X2tX2/S[1,1];
  }
  Out$Ibb = Ibb;

  # Information for alpha
  Iaa = array(0,dim=c(r,r));
  if(n0>0){
    Iaa = Iaa+L[2,2]*P$IPs$Z0tZ0;
  }
  if(n1>0){
    Iaa = Iaa+P$IPs$Z1tZ1/S[2,2];
  }
  Out$Iaa = Iaa;

  # Cross information
  Iba = array(0,dim=c(q,r));
  if(n0>0){
    Iba = L[1,2]*P$IPs$X0tZ0;
  }
  Out$Iba = Iba;
  return(Out);
}

