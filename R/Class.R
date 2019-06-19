#' Bivariate Regression Model
#'
#' @slot Coefficients Regression coefficients.
#' @slot Covariance Outcome covariance matrix.
#' @slot Information Information.
#' @slot Residuals Phenotypic residuals.
#' @name bnr-class
#' @rdname bnr-class
#' @exportClass bnr
setClass(Class="bnr",representation=representation(Coefficients="data.frame",Covariance="matrix",Information="matrix",Residuals="matrix"));

########################
# print Method
########################

#' Print for Bivariate Regression Model
#'
#' @param x A \code{bnr} object.
#' @param ... Unused.
#' @export

print.bnr = function(x,...){
  Coeff = x@Coefficients;
  aux = function(v){
    if(is.numeric(v)){return(signif(v,digits=3))}
    else{return(v)};
  };
  Coeff[] = lapply(Coeff,aux);
  print(Coeff);
};

########################
# show Method
########################

#' Show for Bivariate Regression Model
#' @param object A \code{bnr} object.
#' @rdname bnr-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="bnr"),definition=function(object){print.bnr(x=object)});

########################
# coef Method
########################

#' Extract Coefficients from Bivariate Regression Model
#'
#' @param object A \code{bnr} object.
#' @param ... Unused.
#' @param type Either Target or Surrogate.
#' @export

coef.bnr = function(object,...,type=NULL){
  # Coefficient frame
  Coeff = object@Coefficients;
  if(is.null(type)){
    return(Coeff);
  } else {
    Choices = c("Target","Surrogate");
    if(!(type %in% Choices)){
      stop(paste("Select type from among Target or Surrogate."));
    } else {
      Out = Coeff[Coeff$Outcome==type,];
      return(Out);
    }
  }
};

########################
# resid Method
########################

#' Extract Residuals from Bivariate Regression Model
#'
#' @param object A \code{bnr} object.
#' @param ... Unused.
#' @param type Either Target or Surrogate.
#' @export

residuals.bnr = function(object,...,type=NULL){
  if(is.null(type)){
    return(object@Residuals)
  } else {
    Choices = c("Target","Surrogate");
    if(!(type %in% Choices)){
        stop("Select type from among Target or Surrogate.")
      };
    if(type=="Target"){
      return(object@Residuals[,1]);
    } else {
      return(object@Residuals[,2]);
    }
  }
}

########################
# vcov Method
########################

#' Extract Covariance Matrix from Multivariate Regression Model
#'
#' Returns the either the estimated covariance matrix of the outcome,
#' or the information matrix of the regression coefficients. 
#' 
#' @param object A \code{bnr} object.
#' @param ... Unused.
#' @param type Either Outcome or Information. Default is Information.  
#' @param inv Invert the covariance matrix? Default is FALSE.
#' @export

vcov.bnr = function(object,...,type="Information",inv=F){
  Choices = c("Outcome","Information");
  if(!(type %in% Choices)){stop("Select type from among Information or Outcome.")};
  if(type=="Outcome"){
    Out = object@Covariance;
    if(inv){Out = matInv(Out)};
    return(Out);
  } else {
    Out = object@Information;
    # Invert
    if(inv){Out = matInv(Out);}
    return(Out);
  }
};
