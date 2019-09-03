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

#' @title addSNR_spd
#' @description Adds random noise to an SPD matrix, given a signal to noise ratio.
#' @param num_cov num_cov allows scaling the SNR by the number of covariates (Effective SNR = SNR*num_cov).  
#' @param taper if taper=T, reduces the noise exponentially as distance from the diagonal increases (1/2^d).
#' @export
addSNR_spd <- function(A, SNR, num_cov=1,taper=F) {
  SNR = SNR * num_cov
  V = randsym(sizeR(A,1))
  V = V %*% t(V)
  
  #print(V)
  if(taper) {
    for(row in 1:nrow(A)) {
      for(col in 1:ncol(A)) {
        V[row,col] = V[row,col]/(2^(abs(row-col)+1))
      }
    }
  }
  #print(V)
  
  D = V/norm_TpM_spd(A,V)*(norm_TpM_spd(A,A)/SNR)
  Anew = proj_M_spd(expmap_spd(A,D))
  
  #print(norm_TpM_spd(A,A)/(norm_TpM_spd(A,D)))
  #print(norm_TpM_spd(A,A)/norm_TpM_spd(A,Anew))
  if(!isspd(Anew)) warning("WARNING: random SPD is NOT SPD")
  return(Anew)
}
