# MGLM_Riem
Implementing (as well as augmenting and modifying) the algorithms of Kim et al. 2014 for regressing multiple symmetric positive definite matrices against real valued covariates. The original code from which this repo is based was written in Matlab, and is available from NITRC: [mglm_riem](https://www.nitrc.org/projects/riem_mglm). As well, the surrounding works are hosted on teh University of Wisconsin website: [Hyunwoo J. Kim (2014)](http://pages.cs.wisc.edu/~hwkim/projects/riem-mglm/).

# How to Install

You will need to install from github:

```{r}
library(devtools)
devtools::install_github("mrparker909/MGLMRiem")
```

# How to Use

Work in progress, expect tutorials and examples at a future date. Feel free to contact me through github if you have any questions or concerns!

# Disclaimer
This R package is under development, and bugs can be expected as well as sudden changes to function call formats, function return values, and general package structure. Use at your own risk.

# References

Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352

