test_that("mglm_spd works dim=3x3", {
  X = matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3)
  Y = array(0, dim=c(3,3,5))
  for(i in 1:5) {
    Y[,,i] <- randspd_FAST(3)
  }
  
  mg <- mglm_spd(X = X,Y = Y, maxiter = 50)
  expect_equal(is.null(mg$Yhat), F)
  expect_equal(calc_Rsqr_spd(mg$Y, mg$Yhat)>=0, T)
  expect_equal(calc_Rsqr_spd(mg$Y, mg$Yhat)<=1, T)
})

test_that("mglm_spd works, dim=5x5", {
  set.seed(3)
  d = 4 # dimension of tensor
  X = matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3)
  Y = array(0, dim=c(d,d,5))
  for(i in 1:5) {
    Y[,,i] <- randspd_FAST(d, 15)
  }
  
  mg <- mglm_spd(X = X,Y = Y, maxiter = 50)
  expect_equal(is.null(mg$Yhat), F)
})

test_that("mglm_spd works pKarcher", {
  X = matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3)
  Y = array(0, dim=c(3,3,5))
  for(i in 1:5) {
    Y[,,i] <- randspd_FAST(3)
  }
  
  mg1 <- mglm_spd(X = X,Y = Y, maxiter = 200, pKarcher=T)
  mg2 <- mglm_spd(X = X,Y = Y, maxiter = 200, pKarcher=F)
  
  expect_equal(abs(mg1$E[1]-mg2$E[1]) < 1e-4, T)
})

