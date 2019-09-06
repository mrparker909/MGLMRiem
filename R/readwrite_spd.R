#' @title writeSPD
#' @description write dxdxN array of SPD matrices to disk as a data.frame.
#' @param Y a dxdxN array of SPD matrices
#' @param ID a vector of ID's associated with the SPD matrices (defaults to integers 1:N).
#' @export
writeSPD <- function(Y, ID=NULL, filename="mat.spd") {
  N = sizeR(Y,3)
  if(is.null(ID)) {
    ID = 1:N
  }
  
  siz = nrow(Y[,,1])^2
  Ynrow = nrow(Y[,,1])
  
  df <- data.frame(matrix(ncol = (2+siz), nrow = 0))
  
  # loop over subjects
  for(i in 1:N) {
    subjectID = ID[i]
    Yi = Y[,,i]
    yy = as.numeric(Yi)
    
    df <- rbind(df,data.frame(subjectID,Ynrow,t(yy)))
  }
  write.csv(df, filename, row.names=F)
}

#' @title readSPD
#' @description read a file containing SPD matrices written using writeSPD().
#' @param filename name of file containing SPD matrix data.
#' @export
readSPD <- function(filename="mat.spd") {
  f = read.csv(filename, header=T, row.names=NULL)
  nr = f$Ynrow[1]
  N = nrow(f)
  
  spd = array(0, dim=c(nr,nr,N))
  
  for(i in 1:N) {
    spd[,,i] = as.matrix(f[i,-c(1:(ncol(f)-nr^2))],nrow=nr)
  }
  
  return(list(Y=spd, ID=f$subjectID))
}
