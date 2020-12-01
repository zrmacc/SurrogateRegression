# Surrogate Outcome Regression Analsyis

Zachary McCaw <br>
Updated: 2020-12-01

### Description

This package performs estimation and inference on a partially missing target outcome while borrowing information from a correlated surrogate outcome to increase estimation precision and improve power. The target and surrogate outcomes are jointly modeled within a bivariate outcome regression framework. Unobserved values of either outcome are regarded as missing data. Estimation in the presence of bilateral outcome missingness is performed via an expectation conditional maximization algorithm. A flexible association test is provided for evaluating hypotheses about the target regression parameters. The method and its application to eQTL mapping are described in: [Cross-tissue eQTL mapping in the presence of missing data via surrogate outcome analysis](https://www.biorxiv.org/content/10.1101/2020.11.29.403063v1). Also see:

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
fit <- Fit.BNR(
  t = Y[, "Target"],
  s = Y[, "Surrogate"],
  X = X,
  Z = Z
)
```

```
## Objective increment:  0.0982 
## Objective increment:  0.00172 
## Objective increment:  4.75e-05 
## Objective increment:  1.69e-06 
## Objective increment:  6.54e-08 
## 4 update(s) performed before tolerance limit.
```

```r
show(fit)
```

```
##     Outcome Coefficient  Point     SE      L      U         p
## 1    Target          x1  1.030 0.0359  0.962  1.100 8.55e-182
## 2    Target          x2  0.970 0.0351  0.901  1.040 4.04e-168
## 3 Surrogate          z1 -0.967 0.0334 -1.030 -0.901 6.56e-184
## 4 Surrogate          z2 -0.991 0.0339 -1.060 -0.924 8.90e-188
## 
##         Covariance  Point     SE       L      U
## 1           Target 1.0300 0.0516  0.9360 1.1400
## 2 Target-Surrogate 0.0163 0.0385 -0.0222 0.0548
## 3        Surrogate 1.0100 0.0474  0.9170 1.1000
```
