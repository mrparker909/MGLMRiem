# function d = dist_M_spd(X,Y)
# %DIST_M_SPD returns distance between X and Y on SPD manifold.
# %
# %   d = DIST_M_SPD(X,Y)
# %
# %   See also INNERPROD_TPM_SPD, LOGMAP_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:16:53 $
#   
#   V = logmap_spd(X,Y);
# d = sqrt(innerprod_TpM_spd(V,V,X));
# end

#' @title dis_M_spd
#' @description Calculates the distance on the SPD manifold between two SPD matrices X and Y. \eqn{<Log(X,Y)|Log(X,Y)>_{X}^{0.5}}
#' @param X an SPD matrix
#' @param Y an SPD matrix
#' @export
dist_M_spd <- function(X,Y) {
#DIST_M_SPD returns distance between X and Y on SPD manifold.
#
#   d = DIST_M_SPD(X,Y)
#
#   See also INNERPROD_TPM_SPD, LOGMAP_SPD
#
#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:16:53 $
  
#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $  

   
  V = logmap_spd(X,Y)
  d = sqrt(innerprod_TpM_spd(U = V,V = V,P = X))
  
  return(d)
}

#' @title dis_M_euc
#' @description Calculates the distance in euclidean space between two matrices X and Y. \eqn{d = trace(X'Y)^{0.5}}
#' @param X a matrix
#' @param Y a matrix
#' @export
dist_M_euc <- function(X,Y) {
  #DIST_M_EUC returns distance between X and Y in Euclidean space.
  #
  #   d = DIST_M_EUC(X,Y)
  
  d = sqrt(sum(diag( t(X) %*% Y )))
  
  return(d)
}
