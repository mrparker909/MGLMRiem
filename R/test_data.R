set.seed(12345)
n  <- 4 # number of observations (yi,xi)
m  <- 3 # dim of 3x3 SPD matrices
x1 <- rnorm(n = n, mean =  1, sd = 1) # covariate 1
x2 <- rnorm(n = n, mean = -1, sd = 1) # covariate 2
X  <- matrix(c(x1,x2), byrow = FALSE, nrow = n) # design matrix
 
# generate n SPD(m) matrices (the yi's)
Y <- array(0, dim=c(m,m,n))
for(i in 1:n) {
  Yi <- Matrix::nearPD(matrix(rnorm(x1,abs(x2),n=m*m),nrow=m))$mat
  if(!Matrix::isSymmetric(Yi)) { stop("Yi is not symmetric") }
  Y[,,i] <- as.array(Yi)
}
