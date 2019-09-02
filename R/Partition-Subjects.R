# Purpose: Function to partition subjects by missingness
# Updated: 19/07/31

#' Partition Data by Missingness
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @return List containing the partitioned data.

partSubj = function(t,s,X,Z=NULL){
  if(is.null(Z)){Z=X;}
  # Input structure
  t = matrix(t,ncol=1);
  s = matrix(s,ncol=1);
  # Output structure
  Out = list();

  Out$Dims = list();
  # Dimensions
  q = ncol(X);
  r = ncol(Z);
  p = q+r;

  Out$Dims$p = p;
  Out$Dims$q = q;
  Out$Dims$r = r;

  # Generate indicators
  n = length(t);
  R = array(data=0,dim=c(n,2));
  R[,1] = 1*(!is.na(t));
  R[,2] = 1*(!is.na(s));

  Out$Dims$n = n;
  Out$Inds = list();
  Out$Inds$R = R;

  ## Complete cases
  key = (R[,1]==1)&(R[,2]==1);
  t0 = t[key,,drop=F];
  s0 = s[key,,drop=F];
  X0 = X[key,,drop=F];
  Z0 = Z[key,,drop=F];
  n0 = length(t0);

  # Storage
  Out$Complete = list();
  Out$Complete$t0 = t0;
  Out$Complete$s0 = s0;
  Out$Complete$X0 = X0;
  Out$Complete$Z0 = Z0;
  Out$Dims$n0 = n0;

  # Target missingness
  key = (R[,1]==0)&(R[,2]==1);
  t1 = t[key,,drop=F];
  s1 = s[key,,drop=F];
  X1 = X[key,,drop=F];
  Z1 = Z[key,,drop=F];
  n1 = length(t1);

  # Storage
  Out$TarMiss = list();
  Out$TarMiss$t1 = t1;
  Out$TarMiss$s1 = s1;
  Out$TarMiss$X1 = X1;
  Out$TarMiss$Z1 = Z1;
  Out$Dims$n1 = n1;

  # Surrogate missingness
  key = (R[,1]==1)&(R[,2]==0);
  t2 = t[key,,drop=F];
  s2 = s[key,,drop=F];
  X2 = X[key,,drop=F];
  Z2 = Z[key,,drop=F];
  n2 = length(t2);

  # Storage
  Out$SurMiss = list();
  Out$SurMiss$t2 = t2;
  Out$SurMiss$s2 = s2;
  Out$SurMiss$X2 = X2;
  Out$SurMiss$Z2 = Z2;
  Out$Dims$n2 = n2;

  ## Calculate inner products
  Out$IPs = list();

  # Complete cases
  Out$IPs$X0tX0 = matIP(X0,X0);
  Out$IPs$X0tZ0 = matIP(X0,Z0);
  Out$IPs$Z0tZ0 = matIP(Z0,Z0);
  Out$IPs$X0tT0 = matIP(X0,t0);
  Out$IPs$X0tS0 = matIP(X0,s0);
  Out$IPs$Z0tT0 = matIP(Z0,t0);
  Out$IPs$Z0tS0 = matIP(Z0,s0);

  # Target missing
  Out$IPs$Z1tZ1 = matIP(Z1,Z1);
  Out$IPs$Z1tS1 = matIP(Z1,s1);

  # Surrogate missing
  Out$IPs$X2tX2 = matIP(X2,X2);
  Out$IPs$X2tT2 = matIP(X2,t2);

  # Return
  return(Out);
}
