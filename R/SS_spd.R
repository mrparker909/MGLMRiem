#' @title SST_spd
#' @description Calculate the total sum of squares on the SPD manifold.
#' @param Y     dxdxN array of observed dxd SPD matrices.
#' @param Ybar  Ybar is the Karcher mean of the SPD matrices Y.
#' @param niter number of iterations used to calculate Karcher mean (if Ybar==NULL).
#' @export
SST_spd <- function(Y, Ybar=NULL, niter=200) {
  if(is.null(Ybar)) {
    Ybar = karcher_mean_spd(Y, niter=niter)
  }
  SST = gsqerr_spd(repmat(Ybar,sizeR(Y,3)), Y)
  return(SST)
}


#' @title SSR_spd
#' @description Calculate the regression sum of squares on the SPD manifold.
#' @param Yhat  dxdxN array of estimated dxd SPD matrices.
#' @param Ybar  Ybar is the Karcher mean of the SPD matrices Y.
#' @export
SSR_spd <- function(Yhat, Ybar) {
  SSR = gsqerr_spd(repmat(Ybar,sizeR(Yhat,3)), Yhat)
  return(SSR)
}


#' @title SSE_spd
#' @description Calculate the error sum of squares on the SPD manifold.
#' @param Y     dxdxN array of observed dxd SPD matrices.
#' @param Yhat  dxdxN array of estimated dxd SPD matrices.
#' @export
SSE_spd <- function(Y, Yhat) {
  SSE = gsqerr_spd(Y, Yhat)
  return(SSE)
}