# function V = proj_TpM_spd(V)
# %PROJ_TPM_SPD projects a set of tangent V vectors onto TpM. Symmetrization.
# %
# %   See also MGLM_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 16:59:20 $
#   
#   for i = 1:size(V,3)
#   V(:,:,i) = (V(:,:,i)+V(:,:,i)')/2;
#     end
# end

#' @title proj_TpM_spd
#' @description Projects a set of matrices V onto the tangent space of the SPD manifold (transforms each matrix into a symmetric matrix). 
#' @param V A dxdxN array of dxd matrices to project onto the space of symmetric matrices.
#' @export
proj_TpM_spd <- function(V) {
  # PROJ_TPM_SPD projects a set of tangent V vectors onto TpM. Symmetrization.
  #
  # See also MGLM_SPD
              
  # Hyunwoo J. Kim
  # $Revision: 0.1 $  $Date: 2014/06/23 16:59:20 $

  # Migrated to R by Matthew RP Parker
  # $Revision: 0.2 $  $Date: 2019/06/06 $
  if(length(dim(V))==2) { V <- aug3(V) }
  
  for(i in 1:length(V[1,1,])) {
    V[,,i] <- (V[,,i] + t(V[,,i]))/2
  }
  return(V)
}
