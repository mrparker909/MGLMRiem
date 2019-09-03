# function r = norm_TpM_spd(p,v)
# %NORM_TPM_SPD calculates the norm of tangent vector v in TpM on SPD manifolds.
# %
# %    r  = NORM_TPM_SPD(p,v)
# %
# %   See also DIST_M_SPD, NORM_TPM_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 16:24:38 $
#   
#   r = sqrt(innerprod_TpM_spd(v,v,p));
# end

#' @title norm_TpM_spd
#' @description Calculates the norm of the symmetric matrix v on the SPD manifold at the point p.
#' @param p p is a point on the SPD manifold.
#' @param v v is a tangent vector for which to calculate the norm on the SPD manifold.
#' @export
norm_TpM_spd <- function(p,v) {
# NORM_TPM_SPD calculates the norm of tangent vector v in TpM on SPD manifolds.
#
#    r  = NORM_TPM_SPD(p,v)
#
#   See also DIST_M_SPD, NORM_TPM_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 16:24:38 $
  
#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/07 $ 

  r = sqrt(innerprod_TpM_spd(v,v,p))
  return(r)
}