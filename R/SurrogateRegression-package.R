# Purpose: Package documentation
# Updated: 2022-08-05

#' @useDynLib SurrogateRegression, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' SurrogateRegression: Surrogate Outcome Regression Analysis
#'
#' This package performs estimation and inference on a partially missing target
#' outcome while borrowing information from a correlated surrogate outcome.
#' Rather than regarding the surrogate outcome as a proxy for the target
#' outcome, this package jointly models the target and surrogate outcomes within
#' a bivariate regression framework. Unobserved values of either outcome are
#' treated as missing data. In contrast to imputation-based inference, no
#' assumptions are required regarding the relationship between the target and
#' surrogate outcomes. However, in order for surrogate inference to improve
#' power, the target and surrogate outcomes must be correlated, and the target
#' outcome must be partially missing. The primary estimation function is
#' \code{\link{FitBNR}}. In the case of bilateral missingness, i.e. missingness
#' in both the target and surrogate outcomes, estimation is performed via an
#' expectation conditional maximization either (ECME) algorithm. In the case of
#' unilateral target missingness, estimation is performed using an accelerated
#' least squares procedure. Inference on regression parameters for the target
#' outcome is performed using \code{\link{TestBNR}}.
#'
#' @author Zachary R. McCaw
#' @docType package
#' @name SurrogateRegression
NULL
