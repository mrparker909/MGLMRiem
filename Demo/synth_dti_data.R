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

Yp[,,1] = randspd(3,2,3);
Yp(:,:,2) = randspd(3,2,10);
Yp(:,:,3) = randspd(3,2,10);

% Tangent vectors, geodesic bases.
V = zeros(3,3,npivots);
for j =1:npivots
V(:,:,j) = logmap_spd(Yp(:,:,1),Yp(:,:,j+1));
end

%% Generate Ground Truth Data
Y = zeros(3,3,size(X,2));
for i = 1:size(X,2)
Vtmp = zeros(3,3,1);
for j =1:npivots
Vtmp = Vtmp+ V(:,:,j)*X(j,i);
end
Y(:,:,i) = expmap_spd(Yp(:,:,1),Vtmp);
end

%% Sanity check
notspd = 0 ;
for i=1:size(Y,3)
notspd = notspd + (~isspd(Y(:,:,i)));
end
assert(notspd ==0)
%%
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
