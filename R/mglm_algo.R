#' @title sizeR
#' @description Calculates size of each dimension d of a matrix or array M, returns 1 if d is NA
#' @param M M is a matrix or array object.
#' @param d d is the dimension for which a size is requested (d==NULL returns a vector of all dimsension sizes).
#' @export
sizeR <- function(M, d=NULL) {
  if(is.null(d)) {
    d = dim(M)[1:length(dim(M))]
    return(d)
  }
  
  dims <- dim(M)[d]
  dims[which(is.na(dims))] <- 1
  return(dims)
}

#' @title aug3
#' @description Augment a n x m matrix to a n x m x 1 array, inverse of cont3.
#' @param M M is a matrix which should be augmented to a flat array.
#' @export
aug3 <- function(M) {
  if(length(dim(M))==3) { return(M) }
  if(sizeR(M,3)==1) {
    M <- array(M, dim=c(dim(M)[1:2], 1))  
  }
  return(M)
}

#' @title cont3
#' @description Contract an n x m x 1 array to an n x m matrix, inverse of aug3.
#' @param M M is a flat array which should be contracted down to a matrix.
#' @export
cont3 <- function(M) {
  if(!is.na(dim(M)[3])) {
    return(array(M, dim = dim(M)[1:2]))
  } else {
    return(M)
  }
}

# imported functions
#' @export
Iexpm  <- expm::expm
#' @export
Isqrtm <- pracma::sqrtm
