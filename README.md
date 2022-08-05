---
title: "README"
author: "Zachary R. McCaw"
date: "2022-08-05"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Surrogate Outcome Regression Analsyis

Zachary McCaw <br>
Updated: 2022-08-05

### Description

This package performs estimation and inference on a partially missing target outcome (e.g. gene expression in an inaccessible tissue) while borrowing information from a correlated surrogate outcome (e.g. gene expression in an accessible tissue). Rather than regarding the surrogate outcome as a proxy for the target outcome, the target and surrogate outcomes are jointly modeled within a bivariate regression framework. Unobserved values of either outcome are treated as missing data. In contrast to imputation-based inference, no assumptions are required regarding the relationship between the target and surrogate outcomes. However, in order for surrogate inference to improve efficiency, the target and surrogate outcomes must be correlated, and the target outcome must be partially missing. Estimation in the presence of bilateral outcome missingness is performed via an expectation conditional maximization either ([ECME](https://www.jstor.org/stable/2337067)) algorithm. In the case of unilateral target missingness, estimation is performed using an accelerated least squares procedure. A flexible association test is provided for evaluating hypotheses about the target regression parameters. For additional details, see: McCaw ZR, Gaynor SM, Sun R, Lin X: [Leveraging a surrogate outcome to improve inference on a partially missing target outcome](https://onlinelibrary.wiley.com/doi/10.1111/biom.13629) ([bioRxiv](https://www.biorxiv.org/content/10.1101/2020.11.29.403063v4)). Also see:

* [MGMM](https://github.com/zrmacc/MGMM#missingness-aware-gaussian-mixture-models) for estimation of Gaussian Mixture Models in the presence of missing data. 

### Installation


```r
devtools::install_github(repo = 'zrmacc/SurrogateRegression')
```

### Vignette

Model specification, parameter estimation, and inference on target regression coefficients are discussed in the package vignette, available [here](https://github.com/zrmacc/Spray/blob/master/vignettes/Vignette.pdf).

### Compact Example

In the following, target and surrogate outcome data are generated for 1000 subjects. Each outcome depends on a design matrix with an intercept and a single standard normal covariate. The true target regression coefficient is $\beta = (1, 1)'$ and the true surrogate regression coefficient is $\alpha = (-1, -1)'$; 20% of target outcomes are missing, and 10% of surrogate outcomes are missing.


```r
library(SurrogateRegression)
set.seed(100)

# Data generation.
## Observations.
n <- 1e3

## Target design matrix.
X <- cbind(1, rnorm(n))

## Surrogate design matrix.
Z <- cbind(1, rnorm(n))

## Outcome data.
Y <- rBNR(
  X = X, 
  Z = Z,
  b = c(1, 1), 
  a = c(-1, -1),
  t_miss = 0.2,
  s_miss = 0.1
  )

# Fit bivariate outcome model.
fit <- FitBNR(
  t = Y[, "Target"],
  s = Y[, "Surrogate"],
  X = X,
  Z = Z
)
```

```
## Objective increment:  0.166 
## Objective increment:  0.00497 
## Objective increment:  0.000183 
## Objective increment:  7.15e-06 
## Objective increment:  2.85e-07 
## 4 update(s) performed before tolerance limit.
```

```r
show(fit)
```

```
##     Outcome Coefficient  Point     SE      L      U         p
## 1    Target          x1  0.961 0.0368  0.889  1.030 2.10e-150
## 2    Target          x2  0.952 0.0360  0.882  1.020 1.40e-154
## 3 Surrogate          z1 -1.020 0.0318 -1.080 -0.954 9.53e-225
## 4 Surrogate          z2 -1.100 0.0321 -1.160 -1.040 1.26e-258
## 
##         Covariance    Point     SE       L      U
## 1           Target  1.08000 0.0542  0.9820 1.2000
## 2 Target-Surrogate -0.00478 0.0375 -0.0423 0.0327
## 3        Surrogate  0.90800 0.0428  0.8280 0.9960
```
