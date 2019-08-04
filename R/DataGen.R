# Purpose: Data generation for bivariate normal regression model.
# Updated: 19/07/30

########################
# Data Generation
########################

#' Simulate Bivariate Normal Data with Missingness
#'
#' Function to simulate from a bivariate normal regression model with outcomes
#' missing completely at random.
#'
#' @param X Target design matrix.
#' @param Z Surrogate design matrix.
#' @param b Target regression coefficient.
#' @param a Surrogate regression coefficient.
#' @param mt Target missingness in [0,1].
#' @param ms Surrogate missingness in [0,1].
#' @param S 2x2 target-surrogate covariance matrix.
#'
#' @importFrom mvnfast rmvn
#' @export
#' @return Numeric nx2 matrix. The first column contains the target
#'   outcome, the second contains the surrogate outcome.
#'
#' @examples
#' \dontrun{
#' set.seed(100);
#' # Observations
#' n = 1e3;
#' # Target design
#' X = cbind(1,matrix(rnorm(3*n),nrow=n));
#' # Surrogate design
#' Z = cbind(1,matrix(rnorm(3*n),nrow=n));
#' # Target coefficient
#' b = c(-1,0.1,-0.1,0.1);
#' # Surrogate coefficient
#' a = c(1,-0.1,0.1,-0.1);
#' # Covariance structure
#' S = matrix(c(1,0.5,0.5,1),nrow=2);
#' # Data generation, target and surrogate subject to 10% missingness
#' Y = rBNR(X,Z,b,a,mt=0.1,ms=0.1,S=S);
#' }

rBNR = function(X,Z,b,a,mt=0,ms=0,S){
  # Observations
  n = nrow(X);
  # Linear predictors
  ht = MMP(X,b);
  hs = MMP(Z,a);
  # Residuals
  E = rmvn(n=n,mu=c(0,0),sigma=S);
  # Outcomes
  Y = cbind(ht,hs)+E;

  ## Missingness
  # Target
  Mt = floor(mt*n);
  if(Mt>0){
    Draw = sort(sample(x=n,size=Mt,replace=F));
    Y[Draw,1] = NA;
  }
  # Surrogate
  Ms = floor(ms*n);
  if(Ms>0){
    # Remove subjects with missing target outcome as candidates
    Choices = seq(from=1,to=n)[!is.na(Y[,1])];
    Draw = sort(sample(x=Choices,size=Ms,replace=F));
    Y[Draw,2] = NA;
  }

  ## Output
  colnames(Y) = c("T","S");
  rownames(Y) = seq(1:n);
  return(Y);
}
