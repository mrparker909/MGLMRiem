# function S = weightedsum_mx(mx, w)
# %WEIGHTEDSUM_MX sums matrices with weight w.
# %
# %    S = WEIGHTEDSUM_MX(mx, w)
# %    mx is d x d x N matrices.
# %    w is N weights for matrices. w is a column or row vector.
# %    S is the weighted sum of mx.
# %
# %   See also  MGLM_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 00:13:20 $
#   
#   w = reshape(w,[1 1 length(w)]);
# w = repmat(w, [size(mx,1) size(mx,2) 1]);
# S = sum(mx.*w,3);

#' @export
weightedsum_mx <- function(mx, w) {
  # WEIGHTEDSUM_MX sums matrices with weight w.
  #
  #    S = WEIGHTEDSUM_MX(mx, w)
  #    mx is d x d x N matrices.
  #    w is N weights for matrices. w is a column or row vector.
  #    S is the weighted sum of mx.
  #
  #   See also  MGLM_SPD
  
  #   Hyunwoo J. Kim
  #   $Revision: 0.1 $  $Date: 2014/06/23 00:13:20 $
   
  #   Migrated to R by Matthew RP Parker
  #   $Revision: 0.2 $  $Date: 2019/06/07 $  
  
  #w = reshape(w,[1 1 length(w)]);
  if(dim(mx)[3]<3) {mx <- aug3(mx)}
  
  if(length(w) != dim(mx)[3]) {
    stop("weightedsum_mx(): length(w) != dim(mx)[3]")
  }
  
  print("length(w):")
  print(length(w))
  
  print("dim(mx):")
  print(dim(mx))
  
  for(i in 1:length(w)) {
    mx[,,i] <- w[i] * mx[,,i]
  }
  S <- apply(mx, c(1,2), sum)
  
  #w = repmat(w, [size(mx,1) size(mx,2) 1]);
  #S = sum(mx.*w,3);
  return(S) 
}