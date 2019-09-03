
#' @title isspd_mxstack
#' @description Checks an array of matrices to ensure each matrix is symmetric positive definite.
#' @param Y A dxdxN array of N matrices to check.
#' @export
isspd_mxstack <- function(Y) {
  t = 0
  for(i in 1:sizeR(Y,3)) {
    t = t + (!isspd(Y[,,i]))
  }
  t = !t
  return(t)
}