# Purpose: GLS update for regression parameters
# Updated: 19/08/01

#' Regression Update
#'
#' @param P List of partitioned data. See \code{\link{partSubj}}.
#' @param S Target-surrogate covariance matrix.
#' @export
#' @return List containing the generalized least squares estimates of beta and
#'   alpha.

regUpdate = function(P,S){
  # Dimensions
  p = P$Dims$p;
  q = P$Dims$q;
  r = P$Dims$r;
  n0 = P$Dims$n0;
  n1 = P$Dims$n1;
  n2 = P$Dims$n2;

  # Lambda
  L = matInv(S);

  # Information
  Info = regInfo(P=P,S=S);
  Igg = rbind(cbind(Info$Ibb,Info$Iba),cbind(t(Info$Iba),Info$Iaa));
  Iggi = matInv(Igg);

  # Cross products
  b = array(0,dim=c(q,1));
  a = array(0,dim=c(r,1));
  if(n0>0){
    b = b + L[1,1]*P$IPs$X0tT0 + L[1,2]*P$IPs$X0tS0;
    a = a + L[2,1]*P$IPs$Z0tT0 + L[2,2]*P$IPs$Z0tS0;
  }
  if(n1>0){
    a = a + P$IPs$Z1tS1/S[2,2];
  }
  if(n2>0){
    b = b + P$IPs$X2tT2/S[1,1];
  }

  g = rbind(b,a);

  # GLS Estimate
  gls = as.numeric(MMP(Iggi,g));

  Out = list();
  # Output
  Out$b = gls[1:q];
  Out$a = gls[(q+1):p];
  return(Out);
}

