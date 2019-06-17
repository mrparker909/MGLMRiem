# This is the main algorithm from the conference paper:
# Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352
# Credit goes to the paper authors, blame for errors in code goes to Matthew RP Parker 2019

# returns size of each dimension d, returns 1 if d is NA
#' @export
sizeR <- function(M, d=NULL) {
  if(is.null(d)) {
    d = dim(M)[1:length(dim(M))]
    return(d)
  }
  
  dims <- dim(M)[d]
  dims[which(is.na(dims))] <- 1
  return(dims)
}

# augment a n x m matrix to a n x m x 1 array, inverse of cont3
#' @export
aug3 <- function(M) {
  if(sizeR(M,3)==1) {
    M <- array(M, dim=c(dim(M), 1))  
  }
  return(M)
}

# contract an n x m x 1 array to an n x m matrix, inverse of aug3
#' @export
cont3 <- function(M) {
  if(!is.na(dim(M)[3])) {
    return(array(M, dim = dim(M)[1:2]))
  } else {
    return(M)
  }
}
