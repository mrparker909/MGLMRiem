#' @title Frechet Variance
#' @param Y dxdxn array of dxd spd matrices to find the frechet variance of
#' @param pKarcher if not NULL, the karcher mean from which to center the variance. If NULL, the karcher mean will be calculated.
#' @param maxiter if pKarcher is NULL, maxiter is the maximum number of iterations for calculating the karcher mean.
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

#' @title Calculate R Squared
#' @param Y        The observed spd matrices (an nxnxN array of N spd matrices)
#' @param Yhat     The estimated spd matrices (an nxnxN array of N spd matrices)
#' @param pKarcher The karcher mean of the observed spd matrices, or NULL if the karcher mean is to be calculated
#' @param maxiter  If pKarcher is NULL, maxiter is the maximum number of iterations for calculating the karcher mean.
#' @export
calc_Rsqr_spd <- function(Y, Yhat, pKarcher=NULL, maxiter=200) {
  if(is.null(pKarcher)) {
    pKarcher=karcher_mean_spd(Y,niter=maxiter)
  }
  r2stat_spd(pKarcher,Y,Yhat)
}


#' @title Calculate R Squared
#' @param Y        The observed spd matrices (an nxnxN array of N spd matrices)
#' @param Yhat     The estimated spd matrices (an nxnxN array of N spd matrices)
#' @param pKarcher The karcher mean of the observed spd matrices, or NULL if the karcher mean is to be calculated
#' @param maxiter  If pKarcher is NULL, maxiter is the maximum number of iterations for calculating the karcher mean.
#' @export
calc_Rsqr_euc <- function(Y, Yhat, pKarcher=NULL, maxiter=200) {
  if(is.null(pKarcher)) {
    pKarcher=karcher_mean_spd(Y,niter=maxiter)
  }
  r2stat_euc(pKarcher,Y,Yhat)
}