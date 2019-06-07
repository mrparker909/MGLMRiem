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
  mx <- aug3(mx)
  
  T = array(0, dim=c(sizeR(mx,3),1)) #zeros(size(mx,3),1);
  for(i in 1:sizeR(mx,3)) {
    T[i] = (sum(eigen(mx[,,i])$values <= 0+c ) ==0) && issym(mx[,,i]);
  }
  return(T)
}
    
issym <- function(mx) {
  tol = 0.00001;
  S = array(0, dim=c(sizeR(mx,3), 1)) #zeros(size(mx,3),1);
  
  # augment third dimension
  mx <- aug3(mx)
  
  for(i in 1:sizeR(mx,3)) {
    S[i] = (sum(apply(abs(mx[,,i]-t(mx[,,i])), 2, sum)) < tol)
  }

  T = (sum(S) == sizeR(mx,3))
  return(T)
}
   