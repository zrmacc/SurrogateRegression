---
title: "README"
author: "Zachary R. McCaw"
date: "2019-09-02"
output: 
  html_document: 
    keep_md: TRUE
--- 

# Package Vignette




# Contents

* [Setting](#setting)
* [Example Data](#example-data)
* [Inference](#inference)

# Setting
For each of $n$ independent subjects, suppose two continuous outcomes are potentially observed. Let $t_{i}$ denote the *target* outcome, and let $s_{i}$ denote the *surrogate* outcome. Group the target and surrogate outcomes into a bivariate vector $y_{i} = (t_{i},s_{i})$. For each subject, either the target or the surrogate is potentially missing. However, whether or not an outcome is observed is conditionally independent of its value, given the observed data. Suppose the target mean depends on a vector of covariates $x_{i}$, and the surrogate mean depends on a vector of covariates $z_{i}$:

$$
\mu_{T,i} = E[t_{i}|x_{i}] = x_{i}'\beta \\
\mu_{S,i} = E[s_{i}|z_{i}] = z_{i}'\alpha
$$

Let $\mu_{i}=(\mu_{T,i},\mu_{S,i})$ denote the mean vector. Suppose that $y_{i} = \mu_{i} + \epsilon_{i}$ where $\epsilon_{i}$ is a bivariate normal residual with mean zero and covariance matrix $\Sigma$. This package provides procedures for estimation of model parameters $(\beta,\alpha,\Sigma)$, and inference on the target regression parameters $\beta$. In the case of bilateral (target, surrogate) missingness, estimation is performed via expectation maximization. In the case of unilateral target missingness, estimation is performed via least squares. 

# Example Data
Data are simulated for $n=10^{3}$ subjects. The target `X` and surrogate `Z` design matrices each contain an intercept and three standard normal covariates. The regression coefficient for the target outcome is $\beta = (-1,0.1,-0.1,0)$. The regression coefficient for the surrogate outcome is $\alpha = (1,-0.1,0.1,0)$. The target and surrogate outcome each have unit variance $\Sigma_{TT}=\Sigma_{SS}=1$. The target-surrogate covariance (here, correlation) is $\Sigma_{TS}=\Sigma_{ST}=0.5$. An outcome matrix in which 10% of the target and 20% of the surrogate observations are missing is simulated using `rBNR`. 


```r
library(Spray);
set.seed(100);
# Observations
n = 1e3;
# Target Design
X = cbind(1,matrix(rnorm(3*n),nrow=n));
# Surrogate Design
Z = cbind(1,matrix(rnorm(3*n),nrow=n));
# Target Parameter
b = c(-1,0.1,-0.1,0);
# Surrogate Parameter
a = c(1,-0.1,0.1,0);
# Covariance structure
S = matrix(c(1,0.5,0.5,1),nrow=2);
# Generate data
Y = rBNR(X,Z,b,a,mt=0.1,ms=0.2,S=S);
t = Y[,1];
s = Y[,2];
```

## Formatting Assumptions
The target and surrogate outcome vectors (`t`, `s`) both have length $n$. The unobserved values of the target or surrogate outcome are set to `NA`. The target `X` and surrogate `Z` model matrices are numeric, with all factors and interactions expanded. These matrices contain no missing values. 

<!----------------------------------------------------------------------------------------------------------------------> 

# Inference
Wald and Score tests on $\beta$ are specified using a logical vector `L`. The length of `L` should match the number of columns in the target model matrix `X`. An element of `L` is set to `TRUE` if the regression coefficient for the corresponding column of `X` is zero under $H_{0}$. At least one element of `L` must be `TRUE` (i.e. a test must be specified) and at least one element of `L` must be `FALSE` (i.e. a null model must be estimable). 

Below, various hypothses are tested on the example data. The first is an overall test of $H_{0}:\beta_{1}=\beta_{2}=\beta_{3}=0$, which is false. The second assesses $H_{0}:\beta_{1}=\beta_{2}=0$, which is again false, treating $\beta_{3}$ as a nuisance. The final considers $H_{0}:\beta_{3}=0$, which is true, treating $\beta_{1}$ and $\beta_{2}$ as nuisances. All models include an intercept $\beta_{0}$ under the null.


```r
cat("Joint score test of b1 = b2 = b3 = 0","\n");
signif(test.bnr(t,s,X,Z,L=c(F,T,T,T),report=F,test="Wald"),digits=2);
signif(test.bnr(t,s,X,Z,L=c(F,T,T,T),report=F,test="Score"),digits=2);
cat("\n","Joint score test of b1 = b2 = 0, treating b3 as a nuisance","\n");
signif(test.bnr(t,s,X,Z,L=c(F,T,T,F),report=F,test="Wald"),digits=2);
signif(test.bnr(t,s,X,Z,L=c(F,T,T,F),report=F,test="Score"),digits=2);
cat("\n","Individual score test of b3 = 0, treating b2 and b3 as nuisances","\n");
signif(test.bnr(t,s,X,Z,L=c(F,F,F,T),report=F,test="Wald"),digits=2);
signif(test.bnr(t,s,X,Z,L=c(F,F,F,T),report=F,test="Score"),digits=2);
```

```
## Joint score test of b1 = b2 = b3 = 0 
##    Wald      df       p 
## 2.7e+01 3.0e+00 7.2e-06 
##   Score      df       p 
## 2.6e+01 3.0e+00 1.2e-05 
## 
##  Joint score test of b1 = b2 = 0, treating b3 as a nuisance 
##    Wald      df       p 
## 2.6e+01 2.0e+00 2.0e-06 
##   Score      df       p 
## 2.5e+01 2.0e+00 3.3e-06 
## 
##  Individual score test of b3 = 0, treating b2 and b3 as nuisances 
## Wald   df    p 
## 0.34 1.00 0.56 
## Score    df     p 
##  0.34  1.00  0.56
```

