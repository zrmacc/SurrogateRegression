# Purpose: Master testing function for bivariate normal regression
# Updated: 19/06/18

#' Test Bivariate Normal Regression Model.
#'
#' Performs a test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero in the bivariate normal regression
#' model.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Model matrix.
#' @param Z Surrogate model matrix.
#' @param L Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param test Either Score or Wald. Only Wald is available with LS.
#' @param ... Additional arguments accepted if fitting via EM.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#'
#' @return A numeric vector containing the test statistic, the degrees of
#'   freedom, and a p-value.
#'
#' @examples
#' \dontrun{
#' # See `? rBNR` for data generation
#' # See vignette for test description
#' test.bnr(t=t,s=s,X=X,L=c(F,T,T,T));
#' test.bnr(t=t,s=s,X=X,L=c(F,T,T,F));
#' test.bnr(t=t,s=s,X=X,L=c(F,F,F,T));
#' }

test.bnr = function(t,s,X,Z=NULL,L,test="Wald",...){
  # Input check
  if(!is.vector(t)){stop("A numeric vector is expected for t.")};
  if(!is.vector(s)){stop("A numeric vector is expected for s.")};
  if(!is.matrix(X)){stop("A numeric matrix is expected for X.")};
  if((!is.null(Z))&(!is.matrix(Z))){stop("A numeric matrix is expected for Z.")};
  if((sum(is.na(X))>0)||(sum(is.na(Z)>0))){stop("Missing values are not expected in the covariate matrices.")}
  if(!is.logical(L)){stop("A logical vector is expected for L.")};
  if(!(test%in%c("Score","Wald"))){stop("Please selection either: Score or Wald.")};

  # Determine if s contains missing values, or if Z differs from X.
  EM = (sum(is.na(s))>0)|((!is.null(Z))&(!identical(X,Z)));
  # If missingness occurs in s, apply EM algorithm.
  if(EM){
    if(is.null(Z)){Z=X};
    if(test=="Score"){
      Out = Score.bnem(t=t,s=s,X=X,Z=Z,L=L,...);
    } else {
      Out = Wald.bnem(t=t,s=s,X=X,Z=Z,L=L,...);
    }
  } else {
    # Otherwise, apply the least squares procedure.
    Out = Wald.bnls(t=t,s=s,X=X,L=L);
  }
  # Output
  return(Out);
}
