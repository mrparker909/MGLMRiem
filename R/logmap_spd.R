# function v = logmap_spd(P,X)
# %LOGMAP_SPD maps X on SPD manifold to the tangent space at P.
# %
# %    v = LOGMAP_SPD(P,X)
# %
# %    P, X is a SPD matrix.
# %    v is a symmetric matrix.
# %
# %   See also EXPMAP_SPD, INNERPROD_TPM_SPD, DIST_M_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:20:53 $
#   
#   if norm(P-X) < 1e-18
# v = zeros(size(P));
# return
# end
# 
# [U D] = eig(P);
# g = U*sqrt(D);
# invg = inv(g);
# y = invg*X*invg';
# [V S] = eig(y);
# H = g*V;
# v = H*diag(log(diag(S)))*H';
# 
# %rtX = sqrtm(X);
# %invrtX = inv(rtX);
# %v = rtX*logm(invrtX*Y*invrtX)*rtX;
# % v = (v+v')/2;

#' @title logmap_spd
#' @description Logarithmic Map takes X from the SPD manifold and maps it to the tangent space of symmetric matrices at the point P on the SPD manifold. Analagous to X-P in euclidean space.
#' @param P P is a point on the SPD manifold (a positive definite symmetric matrix).
#' @param X X is an SPD matrix to be mapped to the tangent space at P.
#' @export
logmap_spd <- function(P,X) {
#LOGMAP_SPD maps X on SPD manifold to the tangent space at P.
#
#    v = LOGMAP_SPD(P,X)
#
#    P, X is a SPD matrix.
#    v is a symmetric matrix.
#
#   See also EXPMAP_SPD, INNERPROD_TPM_SPD, DIST_M_SPD
#
#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:20:53 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $  

   
  if(norm(P-X,"2") < 1e-16) { return(array(0, dim=dim(P))) }
  
  if(!isspd(P)) {
    print("P:")
    print(P)
    stop("P was not spd in logmap_spd()")
  }
  
  EIG <- eigen(P,symmetric = T) # eigen(P)
  U   <- EIG$vectors
  D   <- diag(EIG$values)
  
  
  if(any(EIG$values<=0)) {
    D <- ifelse(D < 1e-10, 1e-10, D)
    warning("spd P had non-positive eigenvalues in logmap_spd")
  }
  
  g    = U%*%sqrt(D)
  invg = solve(g)
  y    = invg%*%X%*%t(invg)
  
  EIG <- eigen(y,symmetric = T) # eigen(y)
  V   <- EIG$vectors
  S   <- EIG$values
  
  S = ifelse(S<1e-8, 1e-8, S)
  
  if(any(S < 0 )) { stop("spd y has negative eigenvalues in logmap_spd()") }
  if(any(is.infinite(S))) { warning("spd S has 0 as an eigenvalue in logmap_spd()") }
  
  H = g%*%V
  v = H%*%diag(log(S))%*%t(H)

#rtX = sqrtm(X);
#invrtX = inv(rtX);
#v = rtX*logm(invrtX*Y*invrtX)*rtX;
# v = (v+v')/2;

  return(v)
}
