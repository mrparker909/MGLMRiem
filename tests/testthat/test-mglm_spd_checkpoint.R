test_that("checkpointing works", {
  
  # create test data
  X = 1:10
  Y = array(0, dim=c(2,2,10))
  
  set.seed(123)
  for(i in 1:10) {
    Y[,,i] = randspd_FAST(2)
  }
  
  # create checkpoint file
  res1=mglm_spd(X = as.matrix(t(X)), Y = Y, maxiter = 200, enableCheckpoint = T)
  expect_equal(F, res1$converged)
  
  # load checkpoint file
  load("checkpoint.rda")
  
  # delete checkpoint file
  file.remove("checkpoint.rda")
  
  # run from checkpoint file
  res2=mglm_spd_checkpoint(checkpoint, maxiter = 5, enableCheckpoint=T)
  expect_equal(T, res2$converged)
  
  
  # create checkpoint file
  res1=mglm_spd(X = as.matrix(t(X)), Y = Y, maxiter = 200, enableCheckpoint = T, checkpointPath="./output/")
  expect_equal(F, res1$converged)
  
  # load checkpoint file
  load("./output/checkpoint.rda")
  
  # delete checkpoint file
  file.remove("checkpoint.rda")
  
  # run from checkpoint file
  res2=mglm_spd_checkpoint(checkpoint, maxiter = 5, enableCheckpoint=T, checkpointPath="./output/")
  expect_equal(T, res2$converged)
  
  
})
