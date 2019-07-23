# function Anew = addnoise_spd(A, maxerr)
# V = randsym(size(A,1));
# if norm_TpM_spd(A,V) > maxerr
# V = V/norm_TpM_spd(A,V)*maxerr;
# end
# Anew = expmap_spd(A, V);
# 
# function M = randsym(n)
# M = randn(n);
# M = (M+M')/2;

#' @export
randsym <- function(n) {
  M = matrix(rnorm(n*n), ncol=n) #randn(n)
  M = (M+t(M))/2
  return(M)
}  

#' @export
addnoise_spd <- function(A, maxerr) {
  V = randsym(sizeR(A,1))
  if(norm_TpM_spd(A,V) > maxerr) {
    V = V/norm_TpM_spd(A,V)*maxerr
  }
  Anew = expmap_spd(A, V)
  return(Anew)
}

#' add noise relative to A where maxerr is percentage error (in relation to A)
#' @export
addrelnoise_spd <- function(A, maxerr) {
  V = randsym(sizeR(A,1))
  if(norm_TpM_spd(A,V)/(norm_TpM_spd(A,A)*norm_TpM_spd(A,V)) > maxerr) {
    V = V/norm_TpM_spd(A,V)*(norm_TpM_spd(A,A)*norm_TpM_spd(A,V))*maxerr
  }
  Anew = expmap_spd(A, V)
  return(Anew)
}

#' add noise relative to A where SNR is the signal (A) to noise ratio
#' @export
addSNR_spd <- function(A, SNR) {
  V = randsym(sizeR(A,1))
  V = V/norm_TpM_spd(A,V)*(norm_TpM_spd(A,A)*(runif(1,0,SNR)))
  Anew = expmap_spd(A, V)
  return(Anew)
}