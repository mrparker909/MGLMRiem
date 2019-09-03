# MGLM_Riem
 
Implementing (as well as augmenting and modifying) the algorithms of Kim et al. 2014 for regressing multiple symmetric positive definite matrices against real valued covariates. The original code from which this repo is based was written in Matlab, and is available from NITRC: [mglm_riem](https://www.nitrc.org/projects/riem_mglm). As well, the surrounding works are hosted on the University of Wisconsin website: [Hyunwoo J. Kim (2014)](http://pages.cs.wisc.edu/~hwkim/projects/riem-mglm/).

# How to Install

You will need to install from github:


```r
library(devtools)
devtools::install_github("mrparker909/MGLMRiem")
```

# How to Use

## Generate Some Data

We will use the support package [genDataMGLMRiem](https://github.com/mrparker909/genDataMGLMRiem) to generate some example data to work with.


```r
# devtools::install_github("mrparker909/genDataMGLMRiem")
library(genDataMGLMRiem)

# set seed for reproducibility
set.seed(12345)

# 1) Generate confound data
# Confound 1: 10 observations from a binomial distribution, with probability 0.25, and maximum size 2
C1 = genDataMGLMRiem::SimConfoundData(n = 10, D = rbinom, size = 2, prob = 0.25)

# Confound 2: 10 observations from a standard normal distribution (mean 0, variance 1)
C2 = genDataMGLMRiem::SimConfoundData(n = 10, D = rnorm)

# Put them together to form our X observations:
X = rbind(C1, C2)

# 2) Generate dependent SPD data
# generate 10 3x3 SPD matrices with a covariate effect size between 1 and 2.
simData = genDataMGLMRiem::genSPDdata(N = 10, dims = 3, C = X, minDist = 1, maxDist = 2)
```

## Fit a MGLM Model


```r
library(MGLMRiem)
mod1 = mglm_spd(X = simData$X, Y = simData$Y, pKarcher = T)
```

## Model Diagnostics

We can look at the R squared (variance explained within the SPD manifold):


```r
calc_Rsqr_spd(Y = simData$Y, Yhat = mod1$Yhat)
```

```
## [1] 0.7797107
```

We can also look at the sums of squares:


```r
# Looking at only the upper triangular of the SPD matrices, since they are symmetric

# upper triangular of the mean of the observed SPD matrices
up_kar = extractUT(karcher_mean_spd(simData$Y, niter=100), dims = 3, includeDiagonal = T)

# upper triangular of the estimated SPD matrices
up_hat = extractUT(mod1$Yhat, dims = 3, includeDiagonal = T)

# upper triangular of the observed SPD matrices
up_obs = extractUT(simData$Y, dims = 3, includeDiagonal = T)

# Calculate sum of squared error
SSE = 0
for(i in 1:6) {
  SSE = SSE+sum((up_hat[,i]-up_obs[,i])^2)
}

# Calculate total sum of squares
SST = 0
for(i in 1:6) {
  SST = SST+sum((up_kar[,i]-up_obs[,i])^2)
}

# Calculate regression sum of squares
SSR = 0
for(i in 1:6) {
  SSR = SSR+sum((up_kar[,i]-up_hat[,i])^2)
}
```

So we can easily calculate SSE, SSR, and SST.

And from those we can get at the proportions of explained and unexplained variance:

- Explained Variance = SSR/SST = 0.6301951
- Unexplained Variance = SSE/SST = 0.2674859

# Disclaimer
This R package is under development, and bugs can be expected as well as sudden changes to function call formats, function return values, and general package structure. Use at your own risk. Feel free to contact me through github if you have any questions or concerns!


# References

Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352
