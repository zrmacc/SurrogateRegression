# Purpose: Package documentation
# Updated: 19/07/30

#' @useDynLib Spray, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Spray: Surrogate Phenotype Regression Analysis
#'
#' This package performs estimation and inference for the parameters of a
#' bivariate normal regression model in which elements of the response matrix
#' are missing at random. In the case of bilateral missingness, estimation is
#' performed using \code{\link{fit.bnem}}. In the case of unilaterial
#' missingness, estimation is performed using \code{\link{fit.bnls}}. Inference
#' on regression parameters for the target outcome is performed using
#' \code{\link{test.bnr}}.
#'
#' @author Zachary R. McCaw
#' @docType package
#' @name Spray
NULL
