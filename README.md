# MGLMRiem: Multivariate Generalized Linear Models on Riemann Manifolds


 
Implementing in R (as well as augmenting and modifying) the algorithms of Kim et al. 2014 for regressing multiple symmetric positive definite matrices against real valued covariates. The original code from which this repo is based was written in Matlab, and is available from NITRC: [mglm_riem](https://www.nitrc.org/projects/riem_mglm). As well, the surrounding works are hosted on the University of Wisconsin website: [Hyunwoo J. Kim (2014)](http://pages.cs.wisc.edu/~hwkim/projects/riem-mglm/).

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

We should check that the model converged! If it did not, we would either increase the maximum iterations via maxiters, or implement checkpointing.


```r
# Did the model converge?
mod1$converged
```

```
## [1] TRUE
```

Since the model converged, we can look at the R squared (variance explained within the SPD manifold):


```r
# calculate the Rsquared on the manifold
calc_Rsqr_spd(Y = simData$Y, Yhat = mod1$Yhat)
```

```
## [1] 0.7797107
```

We can also look at the sums of squares:


```r
Ybar = karcher_mean_spd(simData$Y, niter=200)
SSE=SSE_spd(simData$Y, mod1$Yhat)
SST=SST_spd(simData$Y, Ybar)
SSR=SSR_spd(mod1$Yhat, Ybar)
```

So we can easily calculate SSE=9.5953538, SSR= 33.8877806, and SST=43.5579652.

And from those we can get at the proportions of explained and unexplained variance:

- Explained Variance = SSR/SST = 0.7779927
- Unexplained Variance = SSE/SST = 0.2202893
- 'Manifold' Variance = SSM/SST = 1 - SSE/SST - SSR/SST = 0.001718

Note that the 'Manifold' variance is error due to the manifold structure of the response, and the location of the residuals in the tangent space of the manifold (rather than on the manifold itself). Thus it could be considered as part of Explained Variance, and is implicitly included in the R squared calculation: R^2 = (SSR+SSM)/SST = 1-SSE/SST = 0.7797107. 

An alternative, more conservative approach, would be to consider SSM as unexplained, in which case you could define a conservative R^2 as R^2 = SSR/SST = 0.7779927.

# Disclaimer
This R package is under development, and bugs can be expected as well as sudden changes to function call formats, function return values, and general package structure. Use at your own risk. Feel free to contact me through github if you have any questions or concerns!


# References

Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352


