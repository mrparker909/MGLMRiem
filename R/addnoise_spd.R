#' @title randsym
#' @description Creates a random n by n symmetric matrix.
#' @param n dimension of the random symmetric matrix.
#' @export
randsym <- function(n) {
  M = matrix(rnorm(n*n), ncol=n) #randn(n)
  M = (M+t(M))/2
  return(M)
}  


#' @title addSNR_spd
#' @description Adds random noise to an SPD matrix, given a signal to noise ratio. Signal and Noise are measured as matrix content (inner product on the SPD manifold).
#' @param A An SPD matrix which will be noisified.
#' @param SNR Signal to Noise ratio to be used in generating the noise.
#' @param num_cov num_cov allows scaling the SNR by the number of covariates (Effective SNR = SNR*num_cov).  
#' @param taper if taper=T, reduces the noise exponentially as distance from the diagonal increases (1/2^d).
#' @export
addSNR_spd <- function(A, SNR=1, num_cov=1,taper=F) {
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



#' @title addNoise_spd
#' @description Adds random noise to an SPD matrix, given a signal to noise ratio. Signal and Noise are measured as distance on the SPD manifold. If A is the SPD matrix, N is the noise matrix, and I is the identity matrix, then Signal=dist(I,A), Noise=dist(A,N)), and SNR = Signal/Noise.
#' @param A An SPD matrix which will be noisified.
#' @param SNR Signal to Noise ratio to be used in generating the noise.
#' @param returnSNR If T, will return the signal to noise ratio of the original matrix A and the returned matrix
#' @export
addNoise_spd <- function(A, SNR=1, returnSNR=F) {
  In = diag(rep(1,times=sizeR(A,1)))
  dA = dist_M_spd(In, A)
  
  # randomize the signal to noise ratio with mean SNR, std dev sqrt(SNR/)/5
  SNRrnd = 0
  while(SNRrnd < 0.005) {
    SNRrnd = SNR + (rnorm(1,0,sqrt(SNR)/5))
  }
  
  # random symmetric matrix for noise
  N0 = randsym(sizeR(A,1)) # N0 symmetric
  N1 = N0 / dist_M_spd(In,expmap_spd(In,N0))
  N2 = (N1) * dA / (SNRrnd)
  
  N3 = N2
  N = paralleltranslateAtoB_spd(In, A, N3)
  
  Anew = expmap_spd(A,N)
  SNRobs=dist_M_spd(In,A) / dist_M_spd(A,Anew)
  
  if(returnSNR) {
    return(list(A=Anew, SNR=SNRobs))
  } else {
    return(Anew)
  }
}
