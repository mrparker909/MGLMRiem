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
  
  Y_SS = numeric(sizeR(Y,3))
  for(i in 1:sizeR(Y,3)) {
    Y_SS[i] = norm_TpM_spd(Ybar, logmap_spd(Ybar,Y[,,i]))
  }
  
  return(sum(Y_SS^2))
}


#' @title SSR_spd
#' @description Calculate the regression sum of squares on the SPD manifold.
#' @param Yhat  dxdxN array of estimated dxd SPD matrices.
#' @param Ybar  Ybar is the Karcher mean of the SPD matrices Y.
#' @export
SSR_spd <- function(Yhat, Ybar) {
  Y_SS = numeric(sizeR(Yhat,3))
  for(i in 1:sizeR(Yhat,3)) {
    Y_SS[i] = norm_TpM_spd(Ybar, logmap_spd(Ybar,Yhat[,,i]))
  }
  
  return(sum(Y_SS^2))
}


#' @title SSE_spd
#' @description Calculate the error sum of squares on the SPD manifold.
#' @param Y     dxdxN array of observed dxd SPD matrices.
#' @param Yhat  dxdxN array of estimated dxd SPD matrices.
#' @export
SSE_spd <- function(Y, Yhat) {
  Y_SS = numeric(sizeR(Yhat,3))
  for(i in 1:sizeR(Yhat,3)) {
    Y_SS[i] = norm_TpM_spd(Y[,,i], logmap_spd(Y[,,i],Yhat[,,i]))
  }
  
  return(sum(Y_SS^2))
}