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
#' @param SNR Signal to Noise ratio to be used in generating the noise. SNR must be larger 0.
#' @examples 
#' set.seed(623766)
#' A = randspd_FAST(n=5)
#' SNR=.25
#' addNoise_spd(A,SNR=SNR)
#' @export
addNoise_spd <- function(A, SNR=1) {
  d = sizeR(A,1)
  In = diag(rep(1,times=d))
  dA = dist_M_spd(In, A)
  
  Anew = doWhile::doWhile(
    do = {
      attempts = attempts+1
      N0 = randsym(sizeR(A,1)) # N0 symmetric
      N1 = N0 / dist_M_spd(In,expmap_spd(In,N0))
      N2 = (N1) * dA / (SNR)
      N  = paralleltranslateAtoB_spd(In, A, N2)
      Anew = expmap_spd(A,N)
    },
    While = { (innerprod_TpM_spd(Anew,Anew,In) >
                 d*innerprod_TpM_spd(A,A,In)/SNR) & attempts < 1000 },
    Return = {Anew},
    vars = list(A=A, SNR=SNR, d=d, In=In, dA=dA, attempts=0)
  )
  
  if(!isspd(Anew)) { stop("ERROR: in addNoise_spd, new matrix is not SPD")}
  return(Anew)
}
