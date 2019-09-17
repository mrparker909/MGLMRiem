# function xbar = karcher_mean_spd(X, W, niter)
# %KARCHER_MEAN_SPD calculates the intrinsic mean with weight W on SPD manifolds.
# %
# %   xbar = KARCHER_MEAN_SPD(X, [], niter)
# %   xbar = KARCHER_MEAN_SPD(X, W, niter)
# %
# %   W is weights.
# %   X is a set of points on SPD manifolds.
# %   xbar is the Karcher mean of X.
# %   niter is the maximum iterations.
# %
# %   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:20:53 $
#   
#   xbar = X(:,:,1);
# 
# if isempty(W)
# for iter = 1:niter
# phi = mean(logmap_pt2array_spd(xbar,X),3);
# xbar = expmap_spd(xbar, phi);
# if norm(phi) < 1e-18
# break
# end
# end
# else
#   W = W/norm(W,1);
# for iter = 1:niter
# tmp = logmap_pt2array_spd(xbar,X);
# wtmp = zeros(size(tmp));
# for i = 1:size(tmp,3)
# wtmp(:,:,i) = W(i)*tmp(:,:,i);
# end
# phi = sum(wtmp,3);
# xbar = expmap_spd(xbar, phi);
# if norm(phi) < 1e-18
# break
# end
# end
# end

#' @title karcher_mean_spd
#' @description Calculates the Karcher mean with weights W on the SPD manifold.
#' @param X X is a dxdxN array of dxd SPD matrices, for which to calculate the Karcher mean.
#' @param W W is a vector of N weights to apply to the SPD matrices (if W == NULL then equal weights are applied).
#' @param niter niter is the maximum number of iterations in calculating the Karcher mean.
#' @export
karcher_mean_spd <- function(X, W=NULL, niter) {
#KARCHER_MEAN_SPD calculates the intrinsic mean with weight W on SPD manifolds.
#
#   xbar = KARCHER_MEAN_SPD(X, [], niter)
#   xbar = KARCHER_MEAN_SPD(X, W, niter)
#
#   W is weights.
#   X is a set of points on SPD manifolds.
#   xbar is the Karcher mean of X.
#   niter is the maximum iterations.
#
#   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:20:53 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/17 $   
  if(!isspd_mxstack(X)) {
    stop("ERROR: in karcher_mean_spd, X is not an array of SPD matrices")
  }
  xbar = X[,,1]
  
  if(is.null(W)) {
    for(iter in 1:niter) {
      phi = apply(logmap_pt2array_spd(xbar,X), c(1,2), mean)
      xbar = expmap_spd(xbar, phi)
      if(max(svd(phi)$d) < 1e-18) {
        break()
      }
    }
  } else {
    W = W/norm(W,"1")
    for(iter in 1:niter) {
      tmp = logmap_pt2array_spd(xbar,X)
      wtmp = array(0, sizeR(tmp))
      for(i in 1:sizeR(tmp,3)) {
        wtmp[,,i] = W[i]%*%tmp[,,i]
      }
      phi = sum(wtmp,3)
      xbar = expmap_spd(xbar, phi)
      if(norm(phi,"2") < 1e-18) {
        break
      }
    }
  }
  return(xbar)
}