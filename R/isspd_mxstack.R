

#' @export
isspd_mxstack <- function(Y) {
  t = 0
  for(i in 1:sizeR(Y,3)) {
    t = t + (!isspd(Y[,,i]))
  }
  t = !t
  return(t)
}