# DEMO on SPD manifolds

set.seed(222)
print('Start.')
for(i in 1:10) {
  seed = trunc(runif(n = 1, min = 1, 1000))
  
  synth_dti_data
  
  Ybar = karcher_mean_spd(Y,[],500);
  [p, V, E, Yhat, gnorm] = mglm_spd(X,Y);
  r2_iterative  = r2stat_spd(Ybar, Y, Yhat);
  
  [ple, Vle, Ele, Yhatle, logY, Yv ] = mglm_logeuc_spd(X,Y);
  
  r2_logeuc  = r2stat_spd(Ybar, Y, Yhatle);
  fprintf('seed %-4d: r2_iterative %f r2_logeuc %f\n', seed, r2_iterative, r2_logeuc)
}
