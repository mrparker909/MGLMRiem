library(MGLMRiem)
# Parameters
npairs = 10
noise = 1

# Synthesize Ground Truth
X1 = c(seq(0,1,0.25),seq(0,1,0.25))# [0:0.25:1 0:0.25:1];    
X2 = c(0,0,0,0,0,1,1,1,1,1)  #[0 0 0 0 0 1 1 1 1 1 ];
X = rbind(X1, X2)#[X1;X2];

npivots = sizeR(X,1) # Number of points except the base point
Yp = array(0, dim=c(3,3,npivots+1))#zeros(3,3,npivots+1);

Yp[,,1] = randspd(3,2,3)
Yp[,,2] = randspd(3,2,10)
Yp[,,3] = randspd(3,2,10)

# Tangent vectors, geodesic bases.
V = array(0, dim=c(3,3,npivots)) #zeros(3,3,npivots);
for(j in 1:npivots) {
  V[,,j] = logmap_spd(Yp[,,1], Yp[,,j+1])
}

# Generate Ground Truth Data
Y = array(0, dim=c(3,3, sizeR(X,2))) #zeros(3,3,size(X,2));
for(i in 1:sizeR(X,2)) {
  Vtmp = array(0, dim=c(3,3,1))  #zeros(3,3,1)
  for(j in 1:npivots) {
    Vtmp = Vtmp + V[,,j]%*%X[j,i]
  }
  Y[,,i] = expmap_spd(Yp[,,1],Vtmp)
}

# Sanity check
notspd = 0
for(i in 1:sizeR(Y,3)) {
  notspd = notspd + (!isspd(Y[,,i]))
}
if(notspd !=0) stop("generated Y is not an spd!")
#
Xsample = [];
Ysample = zeros(3,3,size(Y,3)*npairs);
isample = 1;
for i = 1:npairs
for j = 1:size(Y,3)
Ysample(:,:,isample) = addnoise_spd(Y(:,:,j),noise);
isample = isample + 1;
end
Xsample =[Xsample X];
end
assert(isspd_mxstack(Ysample) == 1)
X = Xsample;
Y = Ysample;
