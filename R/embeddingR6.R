# function vnew = embeddingR6(p, v) 
# %EMBEDDINGR6 embeds a tangnet vector v in TpM onto R6. This transports v from TpM to T_{I}M.
# %
# %   This will be replaced with the group action for n x n spd matrix. 
# %    
# %   p is a base point. v is a tangent vector.
# %    
# %   vnew = EMBEDDINGR6(p,v)
# %
# %   See also EMBEDDINGR6_VECS, LOGMAP_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:37:53 
# 
# %    step 1
# try
# invp = inv(p);
# catch
# invp = pinv(p); % Numerical problem.
# disp('pinv');
# end
# sqrtinvp= sqrtm(invp);
# %    step 2    
# S = sqrtinvp*v*sqrtinvp;
# %    step 3
# %    v = [Sxx Sxy Sxz Syy Syz Szz]';
# v = symmx2vec(S);
# w = [1  sqrt(2) sqrt(2)   1 sqrt(2) 1]';
# vnew = w.*v;
# end


#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/17 $ 

#' @export
embeddingR6 <- function(p, v) { 
#EMBEDDINGR6 embeds a tangnet vector v in TpM onto R6. This transports v from TpM to T_{I}M.
#
#   This will be replaced with the group action for n x n spd matrix. 
#    
#   p is a base point. v is a tangent vector.
#    
#   vnew = EMBEDDINGR6(p,v)
#
#   See also EMBEDDINGR6_VECS, LOGMAP_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:37:53 

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/17 $ 

  #    step 1
  # invp = tryCatch({
  #   solve(p)
  # }, error = function(e) {
  #   MASS::ginv(p)
  # })
  
  sqrtinvp= Isqrtm(invp)$Binv
  #    step 2    
  S = sqrtinvp%*%v%*%sqrtinvp
  #    step 3
  #    v = [Sxx Sxy Sxz Syy Syz Szz]'
  v = symmx2vec(mx = S)
  w = c(1, sqrt(2), sqrt(2), 1, sqrt(2), 1)
  vnew = w*v

  return(vnew)
}
