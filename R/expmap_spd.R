# function exp_p_x = expmap_spd(P,X)
# %EXPMAP_SPD maps tangent vector X onto SPD manifold.
# %
# %    exp_p_x = EXPMAP_SPD(P,X)
# %
# %    P, exp_p_x is a SPD matrix.
# %    X is a symmetric matrix.
# %    
# %   See also LOGMAP_SPD, KARCHER_MEAN_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:26:13 $
#   
#   if norm(X) < 1e-18
# exp_p_x = P;
# return
# end
# [U D] = eig(P);
# g = U*sqrt(D);
# invg = inv(g);
# Y = invg*X*invg';
#     [V S] = eig(Y);
#     gv = g*V;
#     exp_p_x = gv*diag(exp(diag(S)))*gv';
# 
# %    rtP = sqrtm(P);
# %    invrtP = inv(rtP);
# %    exp_p_v = rtP*expm(invrtP*V*invrtP)*rtP;
# end

#' @export
expmap_spd <- function(P,X) {
#EXPMAP_SPD maps tangent vector X onto SPD manifold.
#
#    exp_p_x = EXPMAP_SPD(P,X)
#
#    P, exp_p_x is a SPD matrix.
#    X is a symmetric matrix.
#    
#   See also LOGMAP_SPD, KARCHER_MEAN_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:26:13 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $ 

  if(max(svd(X)$d) < 1e-18) return(P)

  EIG <- eigen(P)
  U <- EIG$vectors
  D <- diag(EIG$values)

  g <- U*sqrt(D)
  invg = solve(g)
  Y = invg%*%X%*%t(invg)
  
  EIG <- eigen(Y)
  V <- EIG$vectors
  S <- diag(EIG$values)

  gv = g%*%V
  exp_p_x = gv%*%exp(diag(S))%*%t(gv)

#    rtP = sqrtm(P);
#    invrtP = inv(rtP);
#    exp_p_v = rtP*expm(invrtP*V*invrtP)*rtP;

  return(exp_p_x)
}
