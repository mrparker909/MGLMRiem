# function v = symmx2vec(mx)
# %SYMMX2VEC converts matrices mx to vectors v. 
# %
# %   v = symmx2vec(mx)
# %
# %   v is a set of n(n+1)/2 dimensional vectors to n by n matrices.
# %   mx is a set of n x n matrices.
# %
# %   See also INVEMBEDDINGR6, VEC2SYMMX
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 15:09:53 $
#   
#   [ nrow ncol ndata ] = size(mx);
#   v = zeros(nrow*(nrow+1)/2,ndata);
#   k =1;
#   for i=1:ncol
#   for j=i:ncol
#   v(k,:) = squeeze(mx(i,j,:))';
#   k = k + 1;
#   end
#   
#   end
#   end

#' @export
symmx2vec <- function(mx) {
#SYMMX2VEC converts matrices mx to vectors v. 
#
#   v = symmx2vec(mx)
#
#   v is a set of n(n+1)/2 dimensional vectors to n by n matrices.
#   mx is a set of n x n matrices.
#
#   See also INVEMBEDDINGR6, VEC2SYMMX

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 15:09:53 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/17 $ 
  if(is.na(dim(mx)[3])) { mx <- aug3(mx) }
  
  nr    <- dim(mx)[1]
  nc    <- dim(mx)[2]
  ndata <- dim(mx)[3]
  
  v = array(0, dim=c(nr*(nr+1)/2,ndata))
  k =1
  for(i in 1:nc) {
    for(j in i:nc) {
      v[k,] = t(drop(mx[i,j,]))
      k = k + 1
    }
  }
  return(v)
}