# function gsr = gsqerr_spd(X, X_hat)
# %GSQERR_SPD returns the sum of geodesic squared error on SPD manifolds.
# %
# %    gsr = GSQERR_SPD(X, X_hat)
# %
# %
# %   X, X_hat is a set of SPD matrices (dim_X x dim_X x N, where dim_X = size(X,1)).
# %
# %   See also FEVAL_SPD, R2STAT_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $
#   
#   [~, ~, ndata] = size(X);
#   gsr = 0;
#   for idata = 1:ndata
#   gsr = gsr + dist_M_spd(X(:,:,idata),X_hat(:,:,idata))^2;
#   end

#' @export
gsqerr_spd <- function(X, X_hat) {
# GSQERR_SPD returns the sum of geodesic squared error on SPD manifolds.
#
#    gsr = GSQERR_SPD(X, X_hat)
#
#
#   X, X_hat is a set of SPD matrices (dim_X x dim_X x N, where dim_X = size(X,1)).
#
#   See also FEVAL_SPD, R2STAT_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $  

  d <- sizeR(X)
  ndata <- d[3]
  
  gsr = 0
  for(idata in 1:ndata) {
    X     <- aug3(X)
    X_hat <- aug3(X_hat)
    gsr = gsr + dist_M_spd(X[,,idata],X_hat[,,idata])^2
  }
  return(gsr)
}