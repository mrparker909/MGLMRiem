# This is the algorithm from the conference paper:
# Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352
# Credit goes to the paper authors, blame for errors in the code and any additions goes to Matthew RP Parker 2019

#' @title mglm_spd
#' @description This is the MGLM algorithm for performing Riemann Regression on SPD manifolds, developed by: Kim, H. J., Adluru, N., Collins, M. D., Chung, M. K., Bendin, B. B., Johnson, S. C., … Singh, V. (2014). Multivariate General Linear Models (MGLM) on Riemannian Manifolds with Applications to Statistical Analysis of Diffusion Weighted Images. 2014 IEEE Conference on Computer Vision and Pattern Recognition, 2705–2712. https://doi.org/10.1109/CVPR.2014.352. This algorithm was adapted from the publicly available matlab code to R by Matthew RP Parker; blame for errors in the transcription go to Matthew RP Parker 2019. Any issues arising from additions and modifications to the code are also the fault of Matthew RP Parker.
#' @param X       X is a dimX by N matrix of column vectors, each row representing a covariate, each column representing an observation.
#' @param Y       Y is a p by p by N array of p by p symmetric positive definite response matrices.
#' @param maxiter maxiter is the maximum number of iterations before the optimization algorithm stops searching for an optimum. If the algorithm stops before reaching maxiters, then the "converged" variable will be set to TRUE, otherwise it will be set to FALSE.
#' @param pKarcher if TRUE, the estimated base point p will start at the Karcher mean of the observed SPD matrices Y, otherwise the starting point will be the identity matrix.
#' @param enableCheckpoint if TRUE, will create a checkpoint file at the end of each iteration. The checkpoint file may be loaded into R using load(checkpoint.rda), and then mglm_spd_checkpoint(checkpoint) can be run to continue running where the MGLM algorithm left off.
#' @param checkpointPath path to write checkpoint.rda file (if enableCheckpoint=TRUE).
#' @param Memory  Memory is the maximum length of the returned list. This can be useful when the number of iterations needed for convergence is very large. Memory==0 gives 'unlimited' list length.
#' @param tol      tolerance for threshold below which gradient is determined to be zero (stopping condition for algorithm). 
#' @return returns a named list containing the following elements: p (the estimated base point on the manifold), V (the set of estimated covariate coefficient tangent vectors), E (the value of the objective function, which is the sum of squared geodesic error, at each iteration), Yhat (the fitted response values), gnorm (the norm of the gradient at each iteration), converged (a flag indicating whether the algorithm converged before maxiter was reached), MGLMsteps (number of iterations taken by the algorithm).
#' @export
mglm_spd <- function(X, Y, maxiter=500, pKarcher=F, enableCheckpoint=F, checkpointPath="./", Memory=0, tol=1e-10) {
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

# X is nxN (n covariates as rows, N observations as columns)
  
  if(enableCheckpoint) {
    print(paste0("CHECKPOINT ENABLED: setting checkpoint file to ", getwd(),checkpointPath,"checkpoint.rda"))
  }
  
  # if X is a vector of appropriate length, convert X to an array
  if(is.null(dim(X))) {
    if(length(X) == dim(Y)[3]) {
      X = array(X, dim=c(1,length(X)))
    }
  }
  
# ndimX = N = number of covariates
  ndimX = sizeR(X,1) #size(X,1);

# ndimY = mxm
  ndimY = sizeR(Y,1) #size(Y,1);
  ndata = sizeR(X,2) #size(X,2);

  if(ndata != sizeR(Y,3)) { stop('Different number of covariate X and response Y') }


# Initialization
  p <- diag(rep(1, times=ndimY)) #p = Y(:,:,1); 
  if(pKarcher) {
    p = karcher_mean_spd(X=Y,niter=200)
  }
  
  V <- array(data = 0, dim = c(ndimY, ndimY, ndimX)) #V = zeros([ndimY ndimY ndimX]);

# Gradient Descent algorith
# Step size
  c1 = 1

# Safeguard parameter
  c2 = 1

  V <- proj_TpM_spd(V)

  E <- NULL #E = [];
  gnorm <- NULL #gnorm = [];
  # E <- c(E,feval_spd(p,V,X,Y))#E = [E; feval_spd(p,V,X,Y)];
  E <- LimitedMemoryList::LMc(new_el = feval_spd(p,V,X,Y))
  
  step = c1
  for(niter in 1:maxiter) {
    Y_hat = prediction_spd(p,V,X)
    
    if(any(is.na(Y_hat))) {
      stop("element of Y_hat is NA in mglm_spd()")
    }
    if(any(is.null(Y_hat))) {
      stop("element of Y_hat is NULL in mglm_spd()")
    }
    
    J = logmap_vecs_spd(Y_hat, Y) 
    
    if(any(is.na(J))) {
      stop("element of J is NA in mglm_spd()")
    }
    if(any(is.null(J))) {
      stop("element of J is NULL in mglm_spd()")
    }
    
    err_TpM = paralleltranslateAtoB_spd(a = Y_hat, b = p, w = J)
    
    if(any(is.na(err_TpM))) {
      stop("element of err_TpM is NA in mglm_spd()")
    }
    if(any(is.null(err_TpM))) {
      stop("element of err_TpM is NULL in mglm_spd()")
    }
    
    gradp = -1*apply(err_TpM, c(1,2), sum)#-sum(err_TpM,3);

    if(any(is.na(gradp))) {
      stop("element of gradp is NA in mglm_spd()")
    }
    if(any(is.null(gradp))) {
      stop("element of gradp is NULL in mglm_spd()")
    }
    
    # v projection on to tanget space
    gradV = array(0, dim=sizeR(V)) # zeros(size(V));

    # Matrix multiplicaton
    for(iV in 1:sizeR(V,3)) {
      gradV[,,iV] = -weightedsum_mx(err_TpM,X[iV,])
      
      if(any(is.na(gradV[,,iV]))) {
        stop("element of gradV[,,iV] is NA in mglm_spd()")
      }
      if(any(is.null(gradV[,,iV]))) {
        stop("element of gradV[,,iV] is NULL in mglm_spd()")
      }
    }
    

    ns = normVs(p,V = gradV)
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
      V_new = paralleltranslateAtoB_spd(a = p,b = p_new,w = V_new)
      E_new = feval_spd(p_new, V_new, X, Y)

      if(E[1] > E_new) {
        p = p_new
        V = proj_TpM_spd(V_new)
        # E = c(E, E_new)
        E <- LimitedMemoryList::LMc(new_el = E_new, old_list = E, M = ifelse(Memory==0, 1+length(E),Memory))

        if(!is.double(gnorm_new)) {
          stop('Numerical error, gnorm_new is not a double (mglm_spd.R)')
        }

        #gnorm = c(gnorm, gnorm_new)
        gnorm <- LimitedMemoryList::LMc(new_el = gnorm_new, old_list = gnorm, M = ifelse(Memory==0, 1+length(gnorm),Memory))
        
        moved = 1
        step = step*2
        break
      }
    }
    # stopping condition
    if(moved != 1) {
      break 
    } else if(gnorm[1] < tol) {
      break
    #} else if(abs(E[1] - E[2]) < tol[[1]]) {
    #  break
    } else {
      # Checkpoint
      if(enableCheckpoint) {
        checkpoint = list(
          X=X, Y=Y, p=p, V=V, E=E, gnorm=gnorm, niter=niter
        )
        #print(paste0("writing checkpoint to ", getwd(),"/checkpoint.rda"))
        save(checkpoint, file=paste0(checkpointPath,"checkpoint.rda"))
      }
    }
  }

  #E = c(E, feval_spd(p,V,X,Y))
  if(length(E) <= length(gnorm)) {
    E <- LimitedMemoryList::LMc(new_el = feval_spd(p,V,X,Y), old_list = E, M = ifelse(Memory==0, 1+length(E),Memory))
  }
  
  Y_hat = prediction_spd(p,V,X)
  
  #[p, V, E, Y_hat, gnorm]
  return(list(X=X, Y=Y, p=p, V=V, E=E, Yhat=Y_hat, gnorm=gnorm, converged=!(niter>=maxiter), MGLMsteps=niter))
}

# NormVs
#' @export
normVs <- function(p,V) {
  if(length(dim(V))==2) { V = aug3(V) }
  ns <- matrix(0, nrow=sizeR(V,3), ncol=1)
  for(i in 1:sizeR(V,3)) {
    ns[i,1] = norm_TpM_spd(p=p,v=V[,,i])
  }
  if(any(is.na(ns))) {
    stop("element of ns is NA in normVs()")
  }
  if(any(is.null(ns))) {
    stop("element of ns is NULL in normVs()")
  }
  return(ns)
}

# Safeguard
#' @export
safeguard <- function(gradp, gradV, p, c2) {
  ns = normVs(p,gradV)
  normgradv = apply(ns,2,sum) # sum(ns)
  ns = normVs(p,gradp)
  normgradp = apply(ns,2,sum)#sum(ns)
  #norms = [ normgradp normgradv];
  if(is.null(normgradv)) {
    stop("normgradv is NULL in safeguard()")
  }
  if(is.na(normgradv)) {
    stop("normgradv is NA in safeguard()")
  }
  if(is.null(normgradp)) {
    stop("normgradp is NULL in safeguard()")
  }
  if(is.na(normgradp)) {
    stop("normgradp is NA in safeguard()")
  }
  maxnorm = max(normgradp, normgradv)#max(norms);
  if(is.null(c2)) {
    stop("c2 is NULL in safeguard()")
  }
  if(is.na(c2)) {
    stop("c2 is NA in safeguard()")
  }
  if(is.null(maxnorm)) {
    stop("maxnorm is NULL in safeguard()")
  }
  if(is.na(maxnorm)) {
    stop("maxnorm is NA in safeguard()")
  }
  if(maxnorm > c2) {
    gradV = gradV*c2/maxnorm
    gradp = gradp*c2/maxnorm
  }
  return(list(gradp, gradV))
}
