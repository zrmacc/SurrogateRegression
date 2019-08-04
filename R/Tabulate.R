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
  # Convert to log scale
  J = diag(point);
  log.point = log(point);
  log.info = matQF(X=J,A=info);
  # Inverse information
  infoi = matInv(info);
  log.infoi = matInv(log.info);
  # Standard errors
  SE = sqrt(diag(infoi));
  log.SE = sqrt(diag(log.infoi));
  # Critical value
  z = qnorm(p=1-(sig/2));
  # CIs
  L = exp(log.point-z*log.SE);
  U = exp(log.point+z*log.SE);
  # Output
  Out = data.frame("Covariance"=names(point),"Point"=point,"SE"=SE,"L"=L,"U"=U);
  rownames(Out) = seq(1:nrow(Out));
  return(Out);
}
