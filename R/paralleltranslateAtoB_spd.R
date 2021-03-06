# function w_new = paralleltranslateAtoB_spd(a, b, w)
# %PARALLELTRANSLATEATOB_SPD transports a set of tangent vectors w from TaM to
# %TbM.
# %
# %   w_new = PARALLELTRANSLATEATOB_SPD(a, b, w)
# %
# %   a, b are points on SPD matrices. 
# %   w is a set of tangent vectors.
# %   w_new is a set of transported tangent vectors.
# %
# %   See also MGLM_SPD
# 
# %   Hyunwoo J. Kim
# %   $Revision: 0.1 $  $Date: 2014/06/23 17:06:21 $
#   
#   if size(a,3) < size(b,3)
# a = repmat(a, [1 1 size(b,3)]);
# elseif size(a,3) > size(b,3)
# b = repmat(b, [1 1 size(a,3)]);
# end
# 
# if size(b,3) ~= size(w,3)
# % a, b are fixed
# % This changes only w.
# fixab = 1;
# P1 = a;
# P2 = b;
# else
#   fixab = 0;
# end
# 
# w_new = zeros(size(w));
# for i = 1:size(w,3)
# if fixab == 0
# P1 = a(:,:,i);
# P2 = b(:,:,i);
# end
# 
# if norm(P1-P2) < 1e-20
# w_new(:,:,i) = w(:,:,i);
# continue
# end
# 
# %       invP1 = inv(P1);
# %       P12 = sqrtm(invP1*P2*P2*invP1);
# %       T12 = P12\invP1*P2;
# %       B = P1\w(:,:,i);
# %       w_new(:,:,i) = P2*T12'*B*T12;
# w_new(:,:,i) = parallel(P1,P2,w(:,:,i));
# 
# end
# 
# %% symmetrization.
# for i = 1:size(w,3)
# w_new(:,:,i) = (w_new(:,:,i)+w_new(:,:,i)')/2;
# end
# end
# 
# function w_new = parallel(p,q,w)
# rtp = sqrtm(p);
# invrtp = inv(rtp);
# v = logmap_spd(p,q);
# r = expm(invrtp*v/2*invrtp);
# w_new = rtp*r*invrtp*w*invrtp*r*rtp;
# end

# repmat function, repeats a matrix X n times, to make an array that has dimension dim(A) by n
#' @export
repmat = function(X,n) {
  array(X, dim = c(dim(X),n))  
}
  
# repmat function, repeats a matrix X n times, tiled along the columns dimension of X
#' @export
repmat2 = function(X,n) {
  array(rep(X,times=n), dim = c(dim(X)[1], dim(X)[2]*n, dim(X)[-c(1,2)]))  
}

parallel <- function(p,q,w) {
  #if(length(dim(p)) > 2) stop(paste0("dim(p)=",dim(p)))
  
  if(any(is.na(p))) {
    stop("element of p is NA in parallel()")
  }
  if(any(is.null(p))) {
    stop("element of p is NULL in parallel()")
  }
  
  if(any(is.na(q))) {
    stop("element of q is NA in parallel()")
  }
  if(any(is.null(q))) {
    stop("element of q is NULL in parallel()")
  }
  
  if(any(is.na(w))) {
    stop("element of w is NA in parallel()")
  }
  if(any(is.null(w))) {
    stop("element of w is NULL in parallel()")
  }
  
  rtpALL = Isqrtm(p)
  rtp = rtpALL$B #sqrtm(p);
  if(any(is.na(rtp))) {
    stop("element of rtp is NA in parallel()")
  }
  if(any(is.null(rtp))) {
    stop("element of rtp is NULL in parallel()")
  }
  invrtp = rtpALL$Binv  #inv(rtp);
  if(any(is.na(invrtp))) {
    stop("element of invrtp is NA in parallel()")
  }
  if(any(is.null(invrtp))) {
    stop("element of invrtp is NULL in parallel()")
  }
  v = logmap_spd(p,q)
  if(any(is.na(v))) {
    stop("element of v is NA in parallel()")
  }
  if(any(is.null(v))) {
    stop("element of v is NULL in parallel()")
  }
  r = Iexpm(invrtp%*%(v/2)%*%invrtp)
  if(any(r > 10^150)) {
    warning("r in parallel() had very large elements, perhaps Yhat is far from the base point p; this could mean that the base point P is too far from the coefficients V.")
  }
  if(any(is.na(r))) {
    stop("element of r is NA in parallel()")
  }
  if(any(is.null(r))) {
    stop("element of r is NULL in parallel()")
  }
  w_new = rtp%*%r%*%invrtp%*%w%*%invrtp%*%r%*%rtp
  
  if(any(is.na(w_new))) {
    print("invrtp")
    print(invrtp)
    print("w")
    print(w)
    print("r")
    print(r)
    print("rtp")
    print(rtp)
    print("w_new")
    print(w_new)
    stop("element of w_new is NA in parallel()")
  }
  if(any(is.null(w_new))) {
    stop("element of w_new is NULL in parallel()")
  }
  
  return(w_new)
}

#' @title paralleltranslateAtoB_spd
#' @description Performs a parallel transport of the symmetrix matrices w, from point a to point b on the SPD manifold.
#' @param a The starting point on the SPD manifold (a symmetric positive definite matrix).
#' @param b The ending point on the SPD manifold (a symmetric positive definite matrix).
#' @param w The set of tangent vectors (symmetrix matrices) which should be moved via parallel transport from a to b.
#' @export
paralleltranslateAtoB_spd <- function(a, b, w) {
# PARALLELTRANSLATEATOB_SPD transports a set of tangent vectors w from TaM to TbM.
#
#   w_new = PARALLELTRANSLATEATOB_SPD(a, b, w)
#
#   a, b are points on SPD matrices. 
#   w is a set of tangent vectors.
#   w_new is a set of transported tangent vectors.
#
#   See also MGLM_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 17:06:21 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/07 $ 
   a <- aug3(a)
   #b <- aug3(b)
   w <- aug3(w)
   
  if(sizeR(a,3) < sizeR(b,3)) {
    a = repmat(a, sizeR(b,3)) #repmat(a, [1 1 size(b,3)])
  } else if(sizeR(a,3) > sizeR(b,3)) {
    b = repmat(b, sizeR(a,3)) #repmat(b, [1 1 size(a,3)])
  }

  P1 = NULL
  P2 = NULL
  fixab = 0
  if(sizeR(b,3) != sizeR(w,3)) {
    # a, b are fixed
    # This changes only w.
    fixab = 1
    P1 = a;
    P2 = b;
  }

  w_new = array(0,dim = sizeR(w)) #zeros(size(w));
  for(i in 1:sizeR(w,3)) {
    if(fixab == 0) {
      P1 = a[,,i]
      P2 = aug3(b)[,,i]
    }

    P1 <- drop(P1)
    P2 <- drop(P2)
    
    if(max(svd(P1-P2)$d) < 1e-20) { # norm(P1-P2) < 1e-20 {
      w_new[,,i] = w[,,i]
      next()
    }

#       invP1 = inv(P1);
#       P12 = sqrtm(invP1*P2*P2*invP1);
#       T12 = P12\invP1*P2;
#       B = P1\w(:,:,i);
#       w_new(:,:,i) = P2*T12'*B*T12;
    w_new[,,i] = parallel(p=P1,q=P2,w=w[,,i])
  }
  
# symmetrization.
  for(i in 1:sizeR(w,3)) {
    w_new[,,i] = (w_new[,,i]+t(w_new[,,i]))/2
  }
  return(w_new)
}
