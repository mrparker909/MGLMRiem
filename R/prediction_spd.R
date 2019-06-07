# function Yhat = prediction_spd(p,V,X)
# %PREDICTION_SPD predicts phat based on estimate p, V and covariate X.
# %
# %   p is a base point (SPD maxtrix). 
# %   V is a set of tangent vectors (3 x 3 x dimX symmetric matrix).
# %   X is a set of covariates, dimX x N column vectors.
# %   p_hat is the prediction.
# %
# %   See also MGLM_SPD, FEVAL_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 00:13:20 $
#   
#   [ndimX ndata] = size(X);
#   Yhat = zeros(size(p,1),size(p,2),ndata);
#   
#   for i = 1:ndata
#   Vi = zeros(size(p));
#   for j = 1:ndimX
#   Vi = Vi+V(:,:,j)*X(j,i);
#   end
#   Yhat(:,:,i) = expmap_spd(p,Vi);
#   end
#   

prediction_spd <- function(p,V,X) {
#PREDICTION_SPD predicts phat based on estimate p, V and covariate X.
#
#   p is a base point (SPD maxtrix). 
#   V is a set of tangent vectors (3 x 3 x dimX symmetric matrix).
#   X is a set of covariates, dimX x N column vectors.
#   p_hat is the prediction.
#
#   See also MGLM_SPD, FEVAL_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 00:13:20 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $  
  d <- sizeR(X)
  ndimX <- d[1]
  ndata <- d[2]
  
  dp <- sizeR(p)
  Yhat = array(0, dim=c(dp[1],dp[2],ndata))
  
  for(i in 1:ndata) {
    Vi <- array(0,dim=dp);
    Vi <- aug3(Vi)
    for(j in 1:ndimX) {
      Vi <- Vi+V[,,j]*X[j,i]
    }
    Yhat[,,i] = expmap_spd(p,Vi)
  }
  return(Yhat)
}