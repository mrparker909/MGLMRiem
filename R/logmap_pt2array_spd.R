# function V = logmap_pt2array_spd(p,X)
# %LOGMAP_PT2ARRAY_SPD returns logmap(p,Y) for SPD manifolds. This is the
# %faster version of LOGMAP_VECS_SPD for one p not a set of base points P.
# %
# %    V = LOGMAP_PT2ARRAY_SPD(p,X)
# %
# %    p is a SPD matrix. (dim_p x dim_p, where dim_p = size(p,1))
# %    X is a set of SPD matrices (dim_p x dim_p x N).
# %    V is a set of symmetric matrices (dim_p x dim_p x N).
# %
# %   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:50:13 $
#   
#   [U D] = eig(p);
#   g = U*sqrt(D);
#   %invg = inv(g);
#   invg = diag(1./sqrt(diag(D)))*U'; % 1.3 X faster
#   
#   V = zeros(size(X));
#   %% For each data
#   for i = 1:size(X,3)
#   if norm(p-X(:,:,i)) < 1e-18
#   V(:,:,i) = zeros(size(p));
#   continue
#   end
#   y = invg*X(:,:,i)*invg';
#   [U S] = eig(y);
#   H = g*U;
#   V(:,:,i) = H*diag(log(diag(S)))*H';
#   end
# 
logmap_pt2array_spd <- function(p,X) {
#LOGMAP_PT2ARRAY_SPD returns logmap(p,Y) for SPD manifolds. This is the
#faster version of LOGMAP_VECS_SPD for one p not a set of base points P.
#
#    V = LOGMAP_PT2ARRAY_SPD(p,X)
#
#    p is a SPD matrix. (dim_p x dim_p, where dim_p = size(p,1))
#    X is a set of SPD matrices (dim_p x dim_p x N).
#    V is a set of symmetric matrices (dim_p x dim_p x N).
#
#   See also LOGMAP_SPD, LOGMAP_PT2ARRAY_SPD, EXPMAP_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:50:13 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/17 $   
  EIG <- eigen(p,symmetric = T) #eigen(p)
  U   <- EIG$vectors
  D   <- diag(EIG$values)

  g = U%*%sqrt(D)
  invg = diag(1/sqrt(diag(D)))%*%t(U)
  
  V = array(0,sizeR(X))
  # For each data
  for(i in 1:sizeR(X,3)) {
    if(norm(p-X[,,i],"2") < 1e-18) {
      V[,,i] = array(0,sizeR(p))
      next()
    }
    y = invg%*%X[,,i]%*%t(invg)
    EIG <- eigen(y,symmetric = T) #eigen(y)
    U   <- EIG$vectors
    S   <- diag(EIG$values)
    
    # deal with small eigenvalues
    diagS = ifelse(diag(S) < 1e-10, 1e-10, diag(S))
    
    H = g%*%U
    V[,,i] = H%*%diag(log(diagS))%*%t(H)
  }
  return(V)
}

