#' @title frechet_variance_spd
#' @description Calculates the Frechet Variance for a set of SPD matrices.
#' @param Y dxdxn array of dxd SPD matrices for which to find the Frechet Variance.
#' @param pKarcher The Karcher mean from which to center the variance. If NULL, the Karcher mean will be calculated.
#' @param maxiter if pKarcher is NULL, maxiter is the maximum number of iterations for calculating the Karcher mean.
#' @export
frechet_variance_spd <- function(Y, pKarcher=NULL, maxiter=200) {
  if(is.null(pKarcher)) {
    pKarcher=karcher_mean_spd(Y,niter=maxiter)
  }
  fVar = 0
  for(i in 1:sizeR(Y,3)) {
    fVar = fVar + dist_M_spd(pKarcher, Y[,,i])^2
  }
  return(fVar/sizeR(Y,3))
}

#' @title calc_Rsqr_spd
#' @description Calculates R Squared on the SPD manifold.
#' @param Y        The observed SPD matrices (an nxnxN array of N SPD matrices)
#' @param Yhat     The estimated SPD matrices (an nxnxN array of N SPD matrices)
#' @param pKarcher The Karcher mean of the observed SPD matrices, or NULL if the Karcher mean is to be calculated
#' @param maxiter  If pKarcher is NULL, maxiter is the maximum number of iterations for calculating the Karcher mean.
#' @export
calc_Rsqr_spd <- function(Y, Yhat, pKarcher=NULL, maxiter=200) {
  if(is.null(pKarcher)) {
    pKarcher=karcher_mean_spd(Y,niter=maxiter)
  }
  r2stat_spd(pKarcher,Y,Yhat)
}


#' @title calc_Rsqr_euc
#' @description Calculate R Squared in euclidean space.
#' @param Y        The observed matrices (an nxnxN array of N matrices)
#' @param Yhat     The estimated matrices (an nxnxN array of N matrices)
#' @param pKarcher The Karcher mean of the observed spd matrices, or NULL if the Karcher mean is to be calculated.
#' @param maxiter  If pKarcher is NULL, maxiter is the maximum number of iterations for calculating the Karcher mean.
#' @export
calc_Rsqr_euc <- function(Y, Yhat, pKarcher=NULL, maxiter=200) {
  if(is.null(pKarcher)) {
    pKarcher=karcher_mean_spd(Y,niter=maxiter)
  }
  r2stat_euc(pKarcher,Y,Yhat)
}