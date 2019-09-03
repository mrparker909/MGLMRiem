# function T = isspd(mx,varargin)
# %ISSPD check mx is a symmetric positive definite matrix.
# %    This check whether the smallest eigen value is bigger than c.
# %    Default c is epsilon.
# %
# %    Example:
#   %        T = isspd(mx)
#   %        T = isspd(mx,C)
#   %
#   %   See also MGLM_LOGEUC_SPD, MGLM_SPD, PROJ_M_SPD
#   
#   %   Hyunwoo J. Kim
#   %   $Revision: 0.1 $  $Date: 2014/06/23 16:40:17$ 
#     
#     if nargin ==2
#   c = varargin{1};
#   else 
#     c = eps;
#   end
#   % Check matrices are symmetric positive definite.
#   T = zeros(size(mx,3),1);
#   for i=1:size(mx,3)
#   T(i) = (sum(eig(mx(:,:,i)) <= 0+c ) ==0) && issym(mx(:,:,i));
#   end
#   end
#   
#   function [T S] = issym(mx)
#   tol = 0.00001;
#   S = zeros(size(mx,3),1);
#   for i = 1:size(mx,3)
#   S(i) = (sum(sum(abs(mx(:,:,i)-mx(:,:,i)'))) < tol);
#     end
#     T = (sum(S) == size(mx,3));
# end

#' @title isspd
#' @description Checks whether a matrix is symmetric positive definite (SPD) or not. 
#' @param mx Matrix to check.
#' @param c The smallest eigenvalue threshold (eigenvalue larger than c implies non-zero eigenvalue, and SPD matrices have all non-zero eigenvalues).
#' @export
isspd <- function(mx,c=.Machine$double.eps) {
  #ISSPD check mx is a symmetric positive definite matrix.
  #    This check whether the smallest eigen value is bigger than c.
  #    Default c is epsilon.
  #
  #    Example:
  #        T = isspd(mx)
  #        T = isspd(mx,C)
  #
  #   See also MGLM_LOGEUC_SPD, MGLM_SPD, PROJ_M_SPD
  
  #   Hyunwoo J. Kim
  #   $Revision: 0.1 $  $Date: 2014/06/23 16:40:17$ 
     
  #   Migrated to R by Matthew RP Parker
  #   $Revision: 0.2 $  $Date: 2019/06/07 $   
  
  # Check matrices are symmetric positive definite.
  if(length(dim(mx))<3) {mx <- aug3(mx)}
  
  t = array(0, dim=c(sizeR(mx,3),1)) #zeros(size(mx,3),1);
  for(i in 1:sizeR(mx,3)) {
    t[i] = (!any(eigen(mx[,,i],symmetric=T,only.values=T)$values <= c )) && issym(mx[,,i])
    #t[i] = (sum(eigen(mx[,,i])$values <= 0+c ) ==0) && issym(mx[,,i])
  }
  return(all(t == T))
}

#' @title issym
#' @description Checks whether a matrix is symmetric or not. Note that this function takes advantage of the package Rfast if it is installed.
#' @param mx Matrix to check.
#' @export
issym <- function(mx) {
  if(suppressWarnings(require(Rfast, quietly = T))) { return(Rfast::is.symmetric(round(mx, 6))) } else {
    tol = 0.00001
    S = array(0, dim=c(sizeR(mx,3), 1)) #zeros(size(mx,3),1);

    # augment third dimension
    if(length(dim(mx))<3) {mx <- aug3(mx)}

    for(i in 1:sizeR(mx,3)) {
      S[i] = (sum(apply(abs(mx[,,i]-t(mx[,,i])), 2, sum)) < tol)
    }

    T = (sum(S) == sizeR(mx,3))
    return(T)
  }
}    
# issym <- function(mx) {
#   tol = 0.00001;
#   S = array(0, dim=c(sizeR(mx,3), 1)) #zeros(size(mx,3),1);
# 
#   # augment third dimension
#   mx <- aug3(mx)
# 
#   for(i in 1:sizeR(mx,3)) {
#     S[i] = (sum(apply(abs(mx[,,i]-t(mx[,,i])), 2, sum)) < tol)
#   }
# 
#   T = (sum(S) == sizeR(mx,3))
#   return(T)
# }
   