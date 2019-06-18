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
