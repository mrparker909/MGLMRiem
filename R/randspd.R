#' @export
randspd <- function(n, c=3, udist=3) {
# RANDSPD generates n by n random symmatrix positive matrix P.
#
#   P = randspd(n)
#   P = randspd(n,c)
#   P = randspd(n,c,udist)
#
#   c is parameter for variance. Bigger c has bigger variance.
#   udist is the upper bound of distance from I to P w.r.t. GL-invariant
#   measure.
#
#   See also SYNTH_DTI_DATA

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 16:03:38 $

#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $

  P = c*matrix(runif(n*n), ncol=n) # c*(rand(n)-0.5);
  P = P%*%t(P)

  while(dist_M_spd(P,diag(rep(1,n))) > udist) {
    P = c*matrix(runif(n*n), ncol=n) # c*(rand(n)-0.5);
    P = P%*%t(P)
  }
  return(P)
}

#' @title randspd_FAST
#' @description Generate a random Symmetric Positive Definite matrix of dimension n, such that the maximum distance on the manifold is less than maxDist. Makes use of the sparsebnUtils::random.spd function from the library sparsebnUtils.
#' @param n        dimensions of nxn SPD matrix
#' @param maxDist  maximum scale of random component P of SPD matrix (projected from I_nxn with expmap(I,m*P), where m in runif(minDist, maxDist)).
#' @param minDist  minimum scale of random component of P of SPD matrix
#' @param showDist if T, will print the distance of random SPD maxtrix from I_nxn
#' @param NUM      the number of SPD matrices to return (if larger than 1, will return an nxnxNUM array of nxn spd matrices).
#' @param maximumSPDValue maximum value of any entry of the SPD matrix (no limit if NULL)
#' @export
randspd_FAST <- function(n, maxDist=3, showDist=F, NUM=1, minDist=0, maximumSPDValue=NULL) {
  if(minDist==maxDist) {
    warning("minDist == maxDist, setting minDist = 0")
    minDist=0
  }
  # 'central point' on spd cone, identity matrix
  In = diag(rep(1,n))

  Y = array(0, dim=c(n,n,NUM))
  for(i in 1:NUM) {
    while(T) {
      # random SPD
      aDist = runif(1,minDist,maxDist)
      P = expmap_spd(In, sparsebnUtils::random.spd(n)*aDist)
      
      if(!is.null(maximumSPDValue)) {
        while(any(unlist(P) > maximumSPDValue)) {
          aDist = runif(1,minDist,maxDist)
          P = expmap_spd(In, sparsebnUtils::random.spd(n)*aDist)
        }
      }
      
      # control distance from In
      curDist = dist_M_spd(P,In)
      while(curDist > maxDist ) {
        # find the half way point between In and P
        P = karcher_mean_spd(array( c( In , P ) , dim = c(n,n,2)), niter=100)
        #L = logmap_spd(P,In)
        #P = expmap_spd(P,L/2)
        
        curDist = dist_M_spd(P,In)
      }
      if(curDist > minDist) break
    }
    if(showDist) print(curDist)
    if(!isspd(P)) stop("ERROR, randsdp_FAST: GENERATED P IS NOT SPD")
    
    Y[,,i] = P
  }
  if(NUM==1) { Y=Y[,,1] }
  return(Y)
}

