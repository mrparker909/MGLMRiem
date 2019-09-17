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
#' @description Adds random noise to an SPD matrix, given a signal to noise ratio. Signal and Noise are measured as distance on the SPD manifold. If A is the SPD matrix, N is the noise matrix, and I is the identity matrix, then Signal=dist(I,A), Noise=dist(A,N)).
#' @param A An SPD matrix which will be noisified.
#' @param SNR Signal to Noise ratio to be used in generating the noise.
#' @export
addNoise_spd <- function(A, SNR=1) {
  maximumAttempts = 1000
  
  In = diag(rep(1,times=sizeR(A,1)))
  
  m = dist_M_spd(In, A) / SNR
  d = rnorm(1,0,m)
  V1 = randsym(sizeR(A,1))
  V1 = V1 %*% t(V1)
  V2 = V1 / dist_M_spd(In,expmap_spd(In,V1)) * d
  V = expmap_spd(In,V2)
  Vpt = paralleltranslateAtoB_spd(V, A, V2)
  Anew = proj_M_spd(expmap_spd(A,Vpt))
  
  attempts=1
  while(!isspd(Anew) & attempts < maximumAttempts) {
    attempts = attempts+1
    m = dist_M_spd(In, A) / SNR
    d = rnorm(1,0,m)
    V1 = randsym(sizeR(A,1))
    V1 = V1 %*% t(V1)
    V2 = V1 / dist_M_spd(In,expmap_spd(In,V1)) * d
    V = expmap_spd(In,V2)
    Vpt = paralleltranslateAtoB_spd(V, A, V2)
    Anew = proj_M_spd(expmap_spd(A,Vpt))
  }
  
  if(!isspd(Anew)) warning("WARNING: random SPD is NOT SPD")
  return(Anew)
}
