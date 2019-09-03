#' @title mglm_spd_checkpoint
#' @description This is the function for continuing the MGLM algorithm from a checkpoint (useful if the model has not yet converged). 
#' @param checkpoint a checkpoint list object containing: X (the covariate data), Y (the response data), p (the current base point estimate), V (the current coefficient estimates), E (the objective function history), gnorm (the norm of the gradient history), niter (the number of iterations in the checkpoint)
#' @param maxiter maxiter is the maximum number of iterations before the optimization algorithm stops searching for an optimum. If the algorithm stops before reaching maxiters, then the "converged" variable will be set to TRUE, otherwise it will be set to FALSE. Note that maxiter does NOT include the number of iterations loaded from the checkpoint.
#' @param enableCheckpoint if TRUE, will create a checkpoint file at the end of each iteration. The checkpoint file may be loaded into R using load(checkpoint.rda), and then mglm_spd_checkpoint(checkpoint) can be run to continue running MGLM algorithm.
#' @param checkpointPath path to write checkpoint.rda file (if enableCheckpoint=TRUE).
#' @return returns a named list containing the following elements: p (the estimated base point on the manifold), V (the set of estimated covariate coefficient tangent vectors), E (the value of the objective function, which is the sum of squared geodesic error, at each iteration), Yhat (the fitted response values), gnorm (the norm of the gradient at each iteration), converged (a flag indicating whether the algorithm converged before maxiter was reached), MGLMsteps (number of iterations taken by the algorithm).
#' @export
mglm_spd_checkpoint <- function(checkpoint, maxiter=500, enableCheckpoint=T, checkpointPath="./") {
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
  
  
  if(enableCheckpoint) {
    print(paste0("writing checkpoint file to ", getwd(),checkpointPath,"checkpoint.rda"))
  }
  
  X=checkpoint$X
  Y=checkpoint$Y
  p=checkpoint$p
  V=checkpoint$V
  E=checkpoint$E
  gnorm=checkpoint$gnorm
  iter= checkpoint$niter
  maxiter = iter+maxiter
  
  ndimX = sizeR(X,1) #size(X,1);
  
  # ndimY = mxm
  ndimY = sizeR(Y,1) #size(Y,1);
  ndata = sizeR(X,2) #size(X,2);
  
  if(ndata != sizeR(Y,3)) { stop('Different number of covariate X and response Y') }
  
  # Gradient Descent algorith
  # Step size
  c1 = 1
  
  # Safeguard parameter
  c2 = 1
  
  E <- c(E,feval_spd(p,V,X,Y))#E = [E; feval_spd(p,V,X,Y)];
  
  step = c1
  for(niter in iter:maxiter) {
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
    } else {
      # Checkpoint
      if(enableCheckpoint) {
        checkpoint = list(
          X=X, Y=Y, p=p, V=V, E=E, gnorm=gnorm, niter=niter
        )
        save(checkpoint,file=paste0(checkpointPath,"checkpoint.rda"))
      }
    }
  }
  
  E = c(E, feval_spd(p,V,X,Y))
  Y_hat = prediction_spd(p,V,X)
  
  #[p, V, E, Y_hat, gnorm]
  return(list(p=p, V=V, E=E, Yhat=Y_hat, gnorm=gnorm, converged=!(niter>=maxiter), MGLMsteps=niter))
}
