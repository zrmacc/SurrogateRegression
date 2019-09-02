# Purpose: Tabulate regression coefficient estimates
# Updated: 19/08/01

#' Tabulate Regression Coefficients
#'
#' @param point Point estimates.
#' @param info Information matrix.
#' @param sig Significance level.
#'
#' @return Data.table containing the point estimate, standard error, confidence
#'   interval, and Wald p-value.

regTab = function(point,info,sig=0.05){
  # Inverse information
  infoi = matInv(info);
  # Standard errors
  SE = sqrt(diag(infoi));
  # Critical value
  z = qnorm(p=1-(sig/2));
  # CIs
  L = point-z*SE;
  U = point+z*SE;
  # P-values
  p = 2*pnorm(q=abs(point/SE),lower.tail=F);
  # Output
  Out = data.frame("Coefficient"=names(point),"Point"=point,"SE"=SE,"L"=L,"U"=U,"p"=p);
  rownames(Out) = seq(1:nrow(Out));
  return(Out);
}

#' Tabulate Covariance Parameters
#'
#' @param point Point estimates.
#' @param info Information matrix.
#' @param sig Significance level.
#'
#' @return Data.table containing the point estimate, standard error, and
#'   confidence interval.

covTab = function(point,info,sig=0.05){
  # Critical value
  z = qnorm(p=1-(sig/2));
  # Inverse information
  infoi = matInv(info);
  # Standard errors
  SE = sqrt(diag(infoi));
  # CIs for covariance on original scale
  L2 = point[2]-SE[2];
  U2 = point[2]+SE[2];
  
  # CI for variances on log scale
  J = diag(point[c(1,3)]);
  log.point = log(point[c(1,3)]);
  log.info = matQF(X=J,A=info[c(1,3),c(1,3)]);
  log.infoi = matInv(log.info);
  log.SE = sqrt(diag(log.infoi));

  # CIs
  L13 = exp(log.point-z*log.SE);
  U13 = exp(log.point+z*log.SE);
  L = c(L13[1],L2,L13[2]);
  U = c(U13[1],U2,U13[2]);
  # Output
  Out = data.frame("Covariance"=names(point),"Point"=point,"SE"=SE,"L"=L,"U"=U);
  rownames(Out) = seq(1:nrow(Out));
  return(Out);
}
