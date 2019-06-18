test_that("gsqerr_spd works", {
  Y = array(0, dim=c(3,3,5))
  for(i in 1:5) {
    Y[,,i] <- randspd(3)
  }
  
  Y_hat = array(0, dim=c(3,3,5))
  for(i in 1:5) {
    Y_hat[,,i] <- randspd(3)
  }
  
  gsqerr_spd(X = Y, X_hat = Y_hat)
  expect_equal(F, F)
})
