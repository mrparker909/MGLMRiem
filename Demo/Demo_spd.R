# DEMO on SPD manifolds

set.seed(2228)
print('Start.')
for(i in 1:10) {
  seed = trunc(runif(n = 1, min = 1, 1000))
  
  source("./Demo/synth_dti_data.R")
  
  Ybar = karcher_mean_spd(Y,niter = 50) # niter=500
  mglm <- mglm_spd(X,Y, maxiter = 50) # maxiter = default (500)
  p <- mglm$p
  V <- mglm$V
  E <- mglm$E
  Yhat <- mglm$Yhat
  gnorm <- mglm$gnorm
  
  r2_iterative  = r2stat_spd(Y_bar = Ybar, Y = Y, Y_hat = Yhat)
  
  mglm_euc <- mglm_logeuc_spd(X,Y)
  ple = mglm_euc$p
  Vle = mglm_euc$V
  Ele = mglm_euc$E
  Yhatle = mglm_euc$Y_hat
  logY   = mglm_euc$logY
  Yv     = mglm_euc$Yv
  
  r2_logeuc  = r2stat_spd(Y_bar = Ybar, Y = Y, Y_hat = Yhatle)
  print(paste0("seed ", seed,":", "\nr2_iterative = ", r2_iterative, "\nr2_logeuc = ", r2_logeuc))

  
  # checking that Yhat and Yhatle are not too different:
  de <- numeric(100)
  Ydiff <- array(0, dim=c(3,3,100))
  for(i in 1:100) {
    Ydiff[,,i] <- Yhat[,,i] - Yhatle[,,i]
    de[i] <- det(Ydiff[,,i])
  }  

}
