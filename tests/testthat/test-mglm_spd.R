test_that("mglm_spd works dim=3x3", {
  X = matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3)
  Y = array(0, dim=c(3,3,5))
  for(i in 1:5) {
    Y[,,i] <- randspd(3)
  }
  
  mg <- mglm_spd(X = X,Y = Y, maxiter = 50)
  expect_equal(is.null(mg$Yhat), F)
})

test_that("mglm_spd works, dim=5x5", {
  set.seed(3)
  d = 4 # dimension of tensor
  X = matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3)
  Y = array(0, dim=c(d,d,5))
  for(i in 1:5) {
    Y[,,i] <- randspd(d, udist = 15)
  }
  
  mg <- mglm_spd(X = X,Y = Y, maxiter = 50)
  expect_equal(is.null(mg$Yhat), F)
})
