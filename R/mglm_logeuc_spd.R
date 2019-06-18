# function [p, V, E, Y_hat, logY, Yv] = mglm_logeuc_spd(X, Y, varargin)
# %MGLM_LOGEUC_SPD performs MGLM for SPD manifolds by Log Euclidean
# %framework. This also does group action to transport the tangent vector
# %logY to the tangent space at I. Yv is transported logY in T_{I}M. Also
# %this centers X.
# %
# %   [p, V, E, Y_hat] = MGLM_LOGEUC_SPD(X, Y)
# %   [p, V, E, Y_hat, logY, Yv] = MGLM_LOGEUC_SPD(X, Y)
# %   [p, V, E, Y_hat, logY, Yv] =  MGLM_LOGEUC_SPD(X, Y, p, logY, Yv)
# %   [p, V, E, Y_hat, logY, Yv] =  MGLM_LOGEUC_SPD(X, Y, p, logY, Yv, NITER)
# %   has one optional parameter NITER for Karcher mean calculation.  
# %
# %   The result is in p, V, E, Y_hat.
# %
# %   X is dimX x N column vectors
# %   Y is a stack of SPD matrices. 3D arrary 3x3xN.
# %   p is the base point/Karcher mean.
# %   V is a set of tangent vectors (d x d x dimX symmetric matrix).
# %   E is the sum of squared geodesic error.
# %   Y_hat is the prediction.
# %   logY is the Y after logarithm map. logY in TpM.
# %
# %   See also MGLM_SPD, KARCHER_MEAN_SPD, WEIGHTEDSUM_MX,
# %   INVEMBEDDINGR6_VECS, EMBEDDINGR6_VECS, PROJ_M_SPD, GSQERR_SPD
# 
# %   $ Hyunwoo J. Kim $  $ 2016/04/20 15:10:27 (CDT) $
#   %   $ Revision: 0.12 $ 
#   
#   ndimX = size(X,1);
#   ndimY = size(Y,1);
#   ndata =  size(X,2);
#   
#   if ndata ~= size(Y,3)
#   error('Different number of input variables and response variables')
#   end
#   
#   if nargin >= 6
#   niter = varargin{1};
#   else
#     niter = 100;
#   end
#   if nargin <= 2
#   p = karcher_mean_spd(Y, [], niter);
#   end
#   if nargin <= 3
#   logY = logmap_vecs_spd(p, Y);
#   end  
#   % xx xy xz yy yz zz
#   % Center X.
#   Xc = X - repmat(mean(X,2),1,ndata);
#   % Embedding in R6
#   if nargin <= 4
#   Yv = embeddingR6_vecs(p,logY);
#   end
#   L = Yv/Xc;% Peform MGLM
#   logYv_hat = L*Xc; % Prediction.
#   
#   %Transport Y^{\wr} back to TpM
#   V_hat = invembeddingR6_vecs(p,logYv_hat);
#   %Transport V back to TpM
#   V = invembeddingR6_vecs(p,L);
#   
#   Y_hat = zeros(ndimY,ndimY,ndata);
#   for i = 1:ndata
#   %Yhat = expmap(p, Y^{\wr})
#   Y_hat(:,:,i) = expmap_spd(p,V_hat(:,:,i));
#   if ~isspd(Y_hat(:,:,i))
#   disp('Numerical error.');
#   Y_hat(:,:,i) = proj_M_spd(Y_hat(:,:,i)); % For numerical problem.
#   end
#   end
#   E = gsqerr_spd(Y_hat,Y); 
#   end
#   

#' @export
mglm_logeuc_spd <- function(X, Y, p=NULL, logY=NULL, Yv=NULL, niter=100) {
#MGLM_LOGEUC_SPD performs MGLM for SPD manifolds by Log Euclidean
#framework. This also does group action to transport the tangent vector
#logY to the tangent space at I. Yv is transported logY in T_{I}M. Also
#this centers X.
#
#   [p, V, E, Y_hat] = MGLM_LOGEUC_SPD(X, Y)
#   [p, V, E, Y_hat, logY, Yv] = MGLM_LOGEUC_SPD(X, Y)
#   [p, V, E, Y_hat, logY, Yv] =  MGLM_LOGEUC_SPD(X, Y, p, logY, Yv)
#   [p, V, E, Y_hat, logY, Yv] =  MGLM_LOGEUC_SPD(X, Y, p, logY, Yv, NITER)
#   has one optional parameter NITER for Karcher mean calculation.  
#
#   The result is in p, V, E, Y_hat.
#
#   X is dimX x N column vectors
#   Y is a stack of SPD matrices. 3D arrary 3x3xN.
#   p is the base point/Karcher mean.
#   V is a set of tangent vectors (d x d x dimX symmetric matrix).
#   E is the sum of squared geodesic error.
#   Y_hat is the prediction.
#   logY is the Y after logarithm map. logY in TpM.
#
#   See also MGLM_SPD, KARCHER_MEAN_SPD, WEIGHTEDSUM_MX,
#   INVEMBEDDINGR6_VECS, EMBEDDINGR6_VECS, PROJ_M_SPD, GSQERR_SPD

#   $ Hyunwoo J. Kim $  $ 2016/04/20 15:10:27 (CDT) $
#   $ Revision: 0.12 $ 
  
#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/17 $ 
  
  ndimX = sizeR(X,1)
  ndimY = sizeR(Y,1)
  ndata = sizeR(X,2)
  
  if(ndata != sizeR(Y,3)) {
    stop('Different number of input variables and response variables')
  }
  
  if(is.null(p)) {
    p = karcher_mean_spd(Y, niter = niter)
  }
  
  if(is.null(logY)) {
    logY = logmap_vecs_spd(X = p, Y = Y)
  } 
  # xx xy xz yy yz zz
  # Center X.
  Xc = X - repmat2(t(t(apply(X, 1, mean))),ndata)
  # Embedding in R6
  if(is.null(Yv)) {
    Yv = embeddingR6_vecs(p = p,V = logY)
  }
  L = pracma::mrdivide(Yv, Xc) # Yv/Xc # Perform MGLM
  logYv_hat = L%*%Xc # Prediction.
  
  # Transport Y^{\wr} back to TpM
  V_hat = invembeddingR6_vecs(p,V = logYv_hat)
  # Transport V back to TpM
  V = invembeddingR6_vecs(p,L)
  
  Y_hat = array(0, dim=c(ndimY,ndimY,ndata))
  for(i in 1:ndata) {
  #Yhat = expmap(p, Y^{\wr})
    Y_hat[,,i] = expmap_spd(P = p, X = V_hat[,,i])
    if(!isspd(Y_hat[,,i])) {
      print('Numerical error.')
      print(i)
      Y_hat[,,i] = proj_M_spd(Re(Y_hat[,,i])) #proj_M_spd(Y_hat[,,i]) # For numerical problem.
    }
  }
  
  E = gsqerr_spd(X = Y_hat,X_hat = Y)
  
  return(list(p=p, V=V, E=E, Y_hat=Y_hat, logY=logY, Yv=Yv))
}
  