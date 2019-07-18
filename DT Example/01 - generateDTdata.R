# We want to remove confounds from the Diffusion Tensors, 
# using Riemman Manifold Regression

# Y are 3x3 symmetric positive definite matrices (representing Diffusion 
# Tensors for each brain voxel)

# C are the confounds C1, C2, ..., Ck
# Given a set of N Y's (the voxels for a single brain), we want to remove
# the effects of C from the Y's (resulting in 3x3 residual matrices, R, 
# for each SPD matrix Y)

# number of voxels
N <- 100

# number of confounds
k <- 5

# noise to add to ground truth
noise <- 1

# generate N sets of k confounds
C <- array(0, dim=c(k, N))
for(i in 1:k) {
  for(j in 1:N) {
    C[i,j] <- rbinom(n = 1, size = 1, prob = 0.75*i/k)
  }
}

# generate N Diffusion Tensors (no dependence on confounds C)
Y <- array(0, dim=c(3,3,N))
for(i in 1:N) {
  Y[,,i] <- randspd_FAST(n = 3)
}


# generate N Diffusion Tensors (with dependence on confounds C)
npivots = sizeR(C,1) # Number of points except the base point
Yp = array(0, dim=c(3,3,npivots+1))#zeros(3,3,npivots+1);

Yp[,,1] = randspd_FAST(3,2) # base point is Yp[,,1]
for(i in 1:(npivots)) {
  Yp[,,i+1] = randspd_FAST(3,10)
}

# Tangent vectors, geodesic bases.
Vp = array(0, dim=c(3,3,npivots)) #zeros(3,3,npivots);
for(j in 1:npivots) {
  Vp[,,j] = logmap_spd(Yp[,,1], Yp[,,j+1])
}

# Generate Ground Truth Data
Y2 = array(0, dim=c(3,3, sizeR(C,2))) #zeros(3,3,size(X,2));
for(i in 1:sizeR(C,2)) {
  Vtmp = array(0, dim=c(3,3,1))  #zeros(3,3,1)
  for(j in 1:npivots) {
    Vtmp = Vtmp + aug3(Vp[,,j]*C[j,i])
  }
  Y2[,,i] = expmap_spd(P=Yp[,,1],X=Vtmp)
}

# Sanity check
notspd = 0
for(i in 1:sizeR(Y2,3)) {
  notspd = notspd + (!isspd(Y2[,,i]))
}
if(notspd > 0) stop("generated Y2 is not an spd!")
#

# now we need to add some noise
Ysample = array(0, dim=c(3,3, N)) #zeros(3,3,size(Y,3)*npairs);

Xsample = NULL
#Ysample = array(0, dim=c(3,3, sizeR(Y2,3)*N)) #zeros(3,3,size(Y,3)*npairs);
Ysample = array(0, dim=c(3,3, sizeR(Y2,3))) #zeros(3,3,size(Y,3)*npairs);
isample = 1
#for(i in 1:N) {
  for(j in 1:sizeR(Y2,3)) {
    Ysample[,,isample] = addnoise_spd(Y2[,,j],noise)
    isample = isample + 1
  }
  Xsample = cbind(Xsample, C) #[Xsample X];
#}
if(isspd_mxstack(Ysample) != 1) { stop("spd stack contains non-spd matrix") }
Y2 = Ysample
X2 = Xsample








