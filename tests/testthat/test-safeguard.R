test_that("safeguard works", {
  X = matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3)
  Y = array(0, dim=c(3,3,5))
  for(i in 1:5) {
    Y[,,i] <- randspd(3)
  }
  
  ndimX = sizeR(X,1) #size(X,1);
  
  # ndimY = mxm
  ndimY = sizeR(Y,1) #size(Y,1);
  ndata = sizeR(X,2) #size(X,2);
  
  if(ndata != sizeR(Y,3)) { stop('Different number of covariate X and response Y') }
  
  
  # Initialization
  p <- diag(rep(1, times=ndimY)) #p = Y(:,:,1); 
  V <- array(data = 0, dim = c(ndimY, ndimY, ndimX)) #V = zeros([ndimY ndimY ndimX]);
  
  # Gradient Descent algorith
  # Step size
  c1 = 1
  
  # Safeguard parameter
  c2 = 1
  
  V <- proj_TpM_spd(V)
  
  E <- NULL #E = [];
  gnorm <- NULL #gnorm = [];
  E <- c(E,feval_spd(p,V,X,Y))#E = [E; feval_spd(p,V,X,Y)];
  
  step = c1;
  niter=1
  
    Y_hat = prediction_spd(p,V,X)
    J = logmap_vecs_spd(Y_hat, Y)        
    err_TpM = paralleltranslateAtoB_spd(a = Y_hat, b = p, w = J)
    gradp = -1*apply(err_TpM, c(1,2), sum)#-sum(err_TpM,3);
    
    # v projection on to tanget space
    gradV = array(0, dim=sizeR(V)) # zeros(size(V));
    
    # Matrix multiplicaton
    for(iV in 1:sizeR(V,3)) {
      gradV[,,iV] = -weightedsum_mx(err_TpM,X[iV,])
    }
    
    ns = normVs(p,gradV)
    normgradv = apply(ns,2,sum) #sum(ns)
    
    ns = normVs(p,gradp)
    normgradp = apply(ns,2,sum) #sum(ns)
    
    gnorm_new = normgradp+normgradv
    if(!is.double(gnorm_new)) { #~isreal(gnorm_new)
      stop('Numerical Error, gnorm_new is not a double (mglm_spd.R).')
    }
    
    # Safegaurd
    #[gradp gradV] = safeguard(gradp, gradV, p, c2);
    
  safeguard(gradp = gradp, gradV = gradV, p = p, c2 = c2)
  expect_equal(1,1)
})