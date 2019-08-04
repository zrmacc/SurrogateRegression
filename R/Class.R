# Purpose: Class definitions for bivariate normal regression model.
# Updated: 19/07/30

#' Bivariate Regression Model
#'
#' @slot Covariance Residual covariance matrix.
#' @slot Covariance.info Information for covariance parameters.
#' @slot Covariance.tab Table of covariance parameters.
#' @slot Regression.info Information for regression coefficients.
#' @slot Regression.tab Table of regression coefficients.
#' @slot Residuals Outcome residuals.
#' @name bnr-class
#' @rdname bnr-class
#' @exportClass bnr
setClass(Class="bnr",representation=representation(Covariance="matrix",
                                                   Covariance.info="matrix",
                                                   Covariance.tab="data.frame",
                                                   Regression.info="matrix",
                                                   Regression.tab="data.frame",
                                                   Residuals="matrix"));

########################
# print Method
########################

#' Print for Bivariate Regression Model
#'
#' @param x A \code{bnr} object.
#' @param ... Unused.
#' @param type Either Regression or Covariance.
#' @export

print.bnr = function(x,...,type="Regression"){
  aux = function(v){
    if(is.numeric(v)){return(signif(v,digits=3))}
    else{return(v)};
  };
  if(type=="Regression"){
    RegTab = x@Regression.tab;
    RegTab[] = lapply(RegTab,aux);
    print(RegTab);
  } else if(type=="Covariance"){
    CovTab = x@Covariance.tab;
    CovTab[] = lapply(CovTab,aux);
    print(CovTab);
  }
};

########################
# show Method
########################

#' Show for Bivariate Regression Model
#' @param object A \code{bnr} object.
#' @rdname bnr-method
#' @importFrom methods show
setMethod(f="show",signature=c(object="bnr"),
          definition=function(object){
            print.bnr(x=object,type="Regression");
            cat("\n");
            print.bnr(x=object,type="Covariance");
            });

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
  Out = object@Regression.tab;
  if(!is.null(type)){
    Choices = c("Target","Surrogate");
    if(!(type %in% Choices)){stop(paste("Select type from among Target or Surrogate."));};
    Out = Out[Out$Outcome==type,];
  }
  # Return
  return(Out);
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
  Out = NULL;
  if(is.null(type)){
    Out = object@Residuals;
  } else {
    Choices = c("Target","Surrogate");
    if(!(type %in% Choices)){stop("Select type from among Target or Surrogate.")};
    if(type=="Target"){
      Out = object@Residuals[,1];
    } else if(type=="Surrogate"){
      Out = object@Residuals[,2];
    }
  }
  return(Out);
}

########################
# vcov Method
########################

#' Extract Covariance Matrix from Multivariate Regression Model
#'
#' Returns the either the estimated covariance matrix of the outcome, the
#' information matrix for regression coefficients, or the information matrix for
#' covariance parameters.
#'
#' @param object A \code{bnr} object.
#' @param ... Unused.
#' @param type Select "Covariance","Outcome",or "Regression". Default is
#'   "Regression".
#' @param inv Invert the covariance matrix? Default is FALSE.
#' @export

vcov.bnr = function(object,...,type="Regression",inv=F){
  Choices = c("Covariance","Outcome","Regression");
  if(!(type %in% Choices)){stop("Select type from among: Covariance, Outcome, Regression.")};
  Out = NULL;
  # Select matrix to output
  if(type=="Covariance"){
    Out = object@Covariance.info;
  } else if(type=="Outcome"){
    Out = object@Covariance;
  } else if(type=="Regression"){
    Out = object@Regression.info;
  }
  # Inversion
  if(inv){Out = matInv(Out)};
  # Return
  return(Out);
};
