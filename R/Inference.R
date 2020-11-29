# Purpose: Master testing function for bivariate normal regression
# Updated: 2020-11-28.

#' Check Initiation
#' 
#' @param init Optional list of initial parameters for fitting the null model.

CheckInit <- function(init) {
 if (!is.null(init)) {
   if ((!is.list(init)) || 
       is.null(names(init)) || 
       !all(names(init) %in% c("b0", "a0", "sigma0"))
      ) {
     stop("If initial parameter are provided, init should take the form of a list with one or more of these elements: a0, b0, sigma0.")
     }
   }
}

#' Check Test Specification
#' 
#' @param is_zero Logical vector, with as many entires as columns in the target model
#'   matrix, indicating which columns have coefficient zero under the null.
#' @param p Number of columns for the target model matrix.

CheckTestSpec <- function(is_zero, p) {
  
  # Degrees of freedom.
  df0 <- sum(is_zero)
  if (length(is_zero) != p) {
    stop("is_zero should have one entry per column of X.")
  }
  if (df0 == 0) {
    stop("At least 1 entry of is_zero should be TRUE.")
  }
  if (df0 == p) {
    stop("At least 1 entry of is_zero should be FALSE.")
  }
  
}


# -----------------------------------------------------------------------------

#' Test Bivariate Normal Regression Model.
#'
#' Performs a test of the null hypothesis that a subset of the regression
#' parameters for the target outcome are zero in the bivariate normal regression
#' model.
#'
#' @param t Target outcome vector.
#' @param s Surrogate outcome vector.
#' @param X Target model matrix.
#' @param Z Surrogate model matrix.
#' @param is_zero Logical vector, with as many entires as columns in the target
#'   model matrix, indicating which columns have coefficient zero under the
#'   null.
#' @param test Either Score or Wald. Only Wald is available for LS.
#' @param ... Additional arguments accepted if fitting via EM. See
#'   \code{\link{Fit.BNEM}}.
#'
#' @importFrom stats model.matrix pchisq resid vcov
#' @export
#'
#' @return A numeric vector containing the test statistic, the degrees of
#'   freedom, and a p-value.
#'
#' @examples
#' # Generate data.
#' set.seed(100)
#' n <- 1e3
#' X <- cbind(1, rnorm(n))
#' Z <- cbind(1, rnorm(n))
#' data <- rBNR(X = X, Z = Z, b = c(1, 0), a = c(-1, 0), t_miss = 0.1, s_miss = 0.1)
#' 
#' # Test 1st coefficient.
#' wald_test1 <- Test.BNR(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(TRUE, FALSE),
#'   test = "Wald"
#' )
#' 
#' score_test1 <- Test.BNR(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(TRUE, FALSE),
#'   test = "Score"
#' )
#' 
#' # Test 2nd coefficient.
#' wald_test2 <- Test.BNR(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(FALSE, TRUE),
#'   test = "Wald"
#' )
#' 
#' score_test2 <- Test.BNR(
#'   t = data[, 1], 
#'   s = data[, 2], 
#'   X = X, 
#'   Z = Z,
#'   is_zero = c(FALSE, TRUE),
#'   test = "Score"
#' )

Test.BNR <- function(
  t, 
  s, 
  X, 
  Z = NULL, 
  is_zero, 
  test = "Wald",
  ...
) {
  
  # Input checks.
  if ((sum(is.na(X)) > 0) || (sum(is.na(Z) > 0))) {
    stop("Missing values are not expected in the covariate matrices.")
  }
  if (!is.logical(is_zero)) {
    stop("A logical vector is expected for is_zero.")
  }
  if (!(test %in% c("Score", "Wald"))) {
    stop("Please selection either: Score or Wald.")
  }

  # Determine if s contains missing values, or if Z differs from X.
  apply_em <- (sum(is.na(s)) > 0) | ((!is.null(Z)) & (!identical(X, Z)))
  
  # If missingness occurs in s, apply EM algorithm.
  if (apply_em) {
    if (test == "Score") {
      out <- Score.BNEM(t = t, s = s, X = X, Z = Z, is_zero = is_zero, ...)
    } else {
      out <- Wald.BNEM(t = t, s = s, X = X, Z = Z, is_zero = is_zero, ...)
    }
  } else {
    
    # Otherwise, apply the least squares procedure.
    out <- Wald.BNLS(t = t, s = s, X = X, is_zero = is_zero)
  }
  # Output
  return(out)
}
