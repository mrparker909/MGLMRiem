# This is the algorithm from the conference paper:
# Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352
# Credit goes to the paper authors, blame for errors in code goes to Matthew RP Parker 2019

#' @export
mglm_spd <- function(X, Y, maxiter=500) {
# MGLM_SPD performs MGLM on SPD manifolds by interative method.
#
#   [p, V, E, Y_hat, gnorm] = MGLM_SPD(X, Y)
#   [p, V, E, Y_hat, gnorm] = MGLM_SPD(X, Y, MAXITER)
#   has optional parameter MAXITER.  
#
#   The result is in p, V, E, Y_hat.
#
#   X is dimX x N column vectors
#   Y is a stack of SPD matrices. 3D arrary 3x3xN.
#   p is a base point.
#   V is a set of tangent vectors (3 x 3 x dimX symmetric matrix).
#   E is the history of the sum of squared geodesic error.
#   Y_hat is the prediction.
#   gnorm is the history of norm of gradients.
#
#   See also WEIGHTEDSUM_MX, MGLM_LOGEUC_SPD

#   Hyunwoo J. Kim
#   $Revision: 0.1 $  $Date: 2014/06/23 00:13:20 $
#   
#   Migrated to R by Matthew RP Parker
#   $Revision: 0.2 $  $Date: 2019/06/06 $

# X is nxN (n observations as rows, N covariates as columns)
# ndimX = N = number of covariates
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
  for(niter in 1:maxiter) {
    Y_hat = prediction_spd(p,V,X)
    J = logmap_vecs_spd(Y_hat, Y)        
    err_TpM = paralleltranslateAtoB_spd(Y_hat, p, J)
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
    temp  <- safeguard(gradp, gradV, p, c2)
    gradp <- temp[[1]]
    gradV <- temp[[2]]
    
    moved = 0
    for(i in 1:50) {
      step = step*0.5
      # Safegaurd for gradv, gradp
      V_new = V - step*gradV
      p_new = expmap_spd(p,-step*gradp)
      if(!isspd(p_new)) {
        p_new = proj_M_spd(p_new)
      }
      V_new = paralleltranslateAtoB_spd(p,p_new,V_new)
      E_new = feval_spd(p_new, V_new, X, Y)

      if(E[length(E)] > E_new) {
        p = p_new
        V = proj_TpM_spd(V_new)
        E = c(E, E_new)

        if(!is.double(gnorm_new)) {
          stop('Numerical error, gnorm_new is not a double (mglm_spd.R)')
        }

        gnorm = c(gnorm, gnorm_new)
        moved = 1
        step = step*2
        break
      }
    }
    # stopping condition
    if(moved != 1 || gnorm[length(gnorm)] < 1e-10) {
      break 
    }
  }

  E = c(E, feval_spd(p,V,X,Y))
  Y_hat = prediction_spd(p,V,X)
  
  #[p, V, E, Y_hat, gnorm]
  return(list(p=p, E=E, Y_hat=Y_hat, gnorm=gnorm))
}

# NormVs
normVs <- function(p,V) {
  for(i in 1:sizeR(V,3)) {
    ns[i,1] = norm_TpM_spd(p,V[,,i])
  }
}

# Safeguard
safeguard <- function(gradp, gradV, p, c2) {
  ns = normVs(p,gradV)
  normgradv = apply(ns,2,sum) # sum(ns)
  ns = normVs(p,gradp)
  normgradp = apply(ns,2,sum)#sum(ns)
  #norms = [ normgradp normgradv];
  maxnorm = max(normgradp, normgradv)#max(norms);
  if(maxnorm > c2) {
    gradV = gradV*c2/maxnorm
    gradp = gradp*c2/maxnorm
  }
  return(list(gradp, gradV))
}