# Purpose: Package documentation
# Updated: 2020-12-01

#' @useDynLib SurrogateRegression, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' SurrogateRegression: Surrogate Outcome Regression Analysis
#'
#' This package performs estimation and inference on a partially missing target
#' outcome while borrowing information from a correlated surrogate outcome to
#' increase estimation precision and improve power.The primary estimation
#' function is \code{\link{Fit.BNR}}. In the case of bilateral missingness,
#' estimation is performed using \code{\link{Fit.BNEM}}. In the case of
#' unilaterial missingness, estimation is performed using
#' \code{\link{Fit.BNLS}}. Inference on regression parameters for the target
#' outcome is performed using \code{\link{Test.BNR}}.
#'
#' @author Zachary R. McCaw
#' @docType package
#' @name SurrogateRegression
NULL
