# Purpose: Package documentation
# Updated: 19/01/28

#' @useDynLib Spray, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Spray: Surrogate Phenotype Regression Analysis
#'
#' This package performs estimation and inference for the parameters of a
#' bivariate normal regression model in which elements of the response matrix
#' are missing at random. Estimation is performed using \code{\link{fit.bnr}}.
#' Inference on regression parameters for the target outcome is performed using
#' \code{\link{test.bnr}}.
#'
#' @author Zachary R. McCaw
#' @docType package
#' @name Spray
NULL
