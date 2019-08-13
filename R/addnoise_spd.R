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

#' add noise relative to A where SNR is the signal (A) to noise ratio,
#' num_cov allows scaling the SNR by the number of covariates
#' @export
addSNR_spd <- function(A, SNR, num_cov=1) {
  SNR = SNR * num_cov
  V = randsym(sizeR(A,1))
  V = V %*% t(V)
  V = V/norm_TpM_spd(A,V)*(norm_TpM_spd(A,A)*(runif(1,0,1/SNR)))
  Anew = expmap_spd(A,V)
  
  return(Anew)
}