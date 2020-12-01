# Purpose: Class definitions for bivariate normal regression model.
# Updated: 2020-10-25

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
setClass(Class = "bnr", representation = representation(
  Covariance = "matrix",
  Covariance.info = "matrix",
  Covariance.tab = "data.frame",
  Regression.info = "matrix",
  Regression.tab = "data.frame",
  Residuals = "matrix"
))

# -----------------------------------------------------------------------------
# Print Method
# -----------------------------------------------------------------------------

#' Print for Bivariate Regression Model
#'
#' @param x \code{bnr} object.
#' @param ... Unused.
#' @param type Either Regression or Covariance.
#' @export

print.bnr <- function(x, ..., type = "Regression") {
  
  # Function to round numeric columns. 
  aux <- function(v) {
    if (is.numeric(v)) {
      return(signif(v, digits = 3))
    }
    else {
      return(v)
    }
  }
  
  if (type == "Regression") {
    reg_tab <- x@Regression.tab
    reg_tab[] <- lapply(reg_tab, aux)
    print(reg_tab)
  } else if (type == "Covariance") {
    cov_tab <- x@Covariance.tab
    cov_tab[] <- lapply(cov_tab, aux)
    print(cov_tab)
  }
}


# -----------------------------------------------------------------------------
# Show Method.
# -----------------------------------------------------------------------------

#' Show for Bivariate Regression Model
#' 
#' @param object \code{bnr} object.
#' @rdname bnr-method
#' @importFrom methods show
setMethod(
  f = "show", signature = c(object = "bnr"),
  definition = function(object) {
    print.bnr(x = object, type = "Regression")
    cat("\n")
    print.bnr(x = object, type = "Covariance")
  }
)


# -----------------------------------------------------------------------------
# Coef Method
# -----------------------------------------------------------------------------

#' Extract Coefficients from Bivariate Regression Model
#'
#' @param object \code{bnr} object.
#' @param ... Unused.
#' @param type Either Target or Surrogate.
#' @export

coef.bnr <- function(object, ..., type = NULL) {
  out <- object@Regression.tab
  if (!is.null(type)) {
    choices <- c("Target", "Surrogate")
    if (!(type %in% choices)) {
      stop(paste("Select type from among Target or Surrogate."))
    }
    out <- out[out$Outcome == type, ]
  }
  # Return
  return(out)
}


# -----------------------------------------------------------------------------
# Resid Method.
# -----------------------------------------------------------------------------

#' Extract Residuals from Bivariate Regression Model
#'
#' @param object A \code{bnr} object.
#' @param ... Unused.
#' @param type Either Target or Surrogate.
#' @export

residuals.bnr <- function(object, ..., type = NULL) {
  out <- NULL
  if (is.null(type)) {
    out <- object@Residuals
  } else {
    choices <- c("Target", "Surrogate")
    if (!(type %in% choices)) {
      stop("Select type from among Target or Surrogate.")
    }
    if (type == "Target") {
      out <- object@Residuals[, 1]
    } else if (type == "Surrogate") {
      out <- object@Residuals[, 2]
    }
  }
  return(out)
}

# -----------------------------------------------------------------------------
# vcov Method
# -----------------------------------------------------------------------------

#' Extract Covariance Matrix from Bivariate Normal Regression Model
#'
#' Returns the either the estimated covariance matrix of the outcome, the
#' information matrix for regression coefficients, or the information matrix for
#' covariance parameters.
#'
#' @param object \code{bnr} object.
#' @param ... Unused.
#' @param type Select "Covariance","Outcome",or "Regression". Default is
#'   "Regression".
#' @param inv Invert the covariance matrix? Default is FALSE.
#' @export

vcov.bnr <- function(object, ..., type = "Regression", inv = FALSE) {
  choices <- c("Covariance", "Outcome", "Regression")
  if (!(type %in% choices)) {
    stop("Select type from among: Covariance, Outcome, Regression.")
  }
  out <- NULL
  
  # Select covariance matrix.
  if (type == "Covariance") {
    out <- object@Covariance.info
  } else if (type == "Outcome") {
    out <- object@Covariance
  } else if (type == "Regression") {
    out <- object@Regression.info
  }
  
  # Inversion.
  if (inv) {
    out <- matInv(out)
  }
  
  # Output.
  return(out)
}
