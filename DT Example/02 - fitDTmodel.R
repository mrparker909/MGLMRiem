library(MGLMRiem)
set.seed(16424)
source("./DT Example/01 - generateDTdata.R")

# perform MGLM
Ybar = karcher_mean_spd(Y,niter = 500) # niter=500
mglm <- mglm_spd(X = C, Y = Y, maxiter = 500) # maxiter = default (500)
p <- mglm$p
V <- mglm$V
E <- mglm$E
Yhat <- mglm$Yhat
gnorm <- mglm$gnorm

# expect a low r2, since Y is independent of X
r2_iterative  = r2stat_spd(Y_bar = Ybar, Y = Y, Y_hat = Yhat)
# r2 = 0.0406

# perform MGLM
Ybar2 = karcher_mean_spd(Y2,niter = 500) # niter=500
mglm2 <- mglm_spd(X = X2, Y = Y2, maxiter = 500) # maxiter = default (500)
p2 <- mglm2$p
V2 <- mglm2$V
E2 <- mglm2$E
Yhat2  <- mglm2$Yhat
gnorm2 <- mglm2$gnorm

# expect a high r2, since Y is dependent on X
r2_iterative2 = r2stat_spd(Y_bar = Ybar2, Y = Y2, Y_hat = Yhat2)
# r2=0.9973