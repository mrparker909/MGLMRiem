# function p = proj_M_spd(X,varargin)
# %PROJ_M_SPD projects a matrix onto SPD manifolds.
# %
# %
# %    Example:
#   %        p = PROJ_M_SPD(X)
#   %
#   %   p is the point on SPD manifolds.
#   %   X is a n x n matrix.
#   %
#   %   See also MGLM_LOGEUC_SPD, MGLM_SPD
#   
#   %   Hyunwoo J. Kim
#   %   $Revision: 0.1 $  $Date: 2014/06/23 16:40:17$ 
#     
#     if nargin == 2
#   c = varargin{1};
#   else 
#     c = eps;
#   end
#   
#   % Make a matrix symmetric positive definite.
#   if norm(X-X') > eps
#           X = (X+X')/2;
#   end
#   [V D ] = eig(X);
#   D = diag(D);
#   p = zeros(size(X));
#   for i =1:length(D)
#   if D(i) > 0
#   p = p + D(i)*V(:,i)*V(:,i)';
#   end
#   end
#   % Now X is spd
#   % Make psd matrix
#   if sum(D > 0+c) < length(D)
#   a = 1e-16; 
#   pnew = p;
#   while ~isspd(pnew, c)
#   pnew = p + a*eye(3);
#   a = 2*a;
#   end
#   p = pnew;
#   end
#   end

#' @export
proj_M_spd <- function(X,c=.Machine$double.eps) {
# PROJ_M_SPD projects a matrix onto SPD manifolds.
#
#
#    Example:
  #        p = PROJ_M_SPD(X)
  #
  #   p is the point on SPD manifolds.
  #   X is a n x n matrix.
  #
  #   See also MGLM_LOGEUC_SPD, MGLM_SPD
  
  #   Hyunwoo J. Kim
  #   $Revision: 0.1 $  $Date: 2014/06/23 16:40:17$ 
    
  #   Migrated to R by Matthew RP Parker
  #   $Revision: 0.2 $  $Date: 2019/06/07 $ 
  
  # Make a matrix symmetric positive definite.
  if(max(svd(X-t(X))$d) > c) {
          X = (X+t(X))/2
  }
  #[V D ] = eig(X);
  EIG <- suppressWarnings(RSpectra::eigs_sym(X,k=dim(X)[1])) # eigen(X)
  V   <- EIG$vectors
  D   <- EIG$values

  p = array(0, dim=sizeR(X)) #zeros(size(X));
  for(i in 1:length(D)) {
    if(D[i] > 0) {
      p = p + D[i]*V[,i]%*%t(V[,i])
    }
  }
  # Now X is spd
  # Make psd matrix
  if(sum(D > c) < length(D)) {
    a = 1e-16 
    pnew = p
    while(!isspd(pnew, c)) {
      pnew = p + a*diag(rep(1,dim(p)[1]))
      a = 2*a
    }
    p = pnew
  }

  return(p)
}