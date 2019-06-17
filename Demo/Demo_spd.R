# DEMO on SPD manifolds

set.seed(222)
print('Start.')
for(i in 1:10) {
  seed = trunc(runif(n = 1, min = 1, 1000))
  
  source("./Demo/synth_dti_data.R")
  
  # continue converting to R from here:
  Ybar = karcher_mean_spd(Y,niter = 500)
  mglm <- mglm_spd(X,Y)
  p <- mglm$p
  V <- mglm$V
  E <- mglm$E
  Yhat <- mglm$Yhat
  gnorm <- mglm$gnorm
  
  r2_iterative  = r2stat_spd(Ybar, Y, Yhat)
  mglm_euc <- mglm_logeuc_spd(X,Y)
  ple = mglm_euc$p
  Vle = mglm_euc$V
  Ele = mglm_euc$E
  Yhatle = mglm_euc$Yhat
  logY   = mglm_euc$logY
  Yv     = mglm_euc$Yv
  
  r2_logeuc  = r2stat_spd(Ybar, Y, Yhatle)
  print(paste0("seed ", seed,":", "\nr2_iterative = ", r2_iterative, "\nr2_logeuc = ", r2_logeuc))
}
