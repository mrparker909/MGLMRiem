test_that("karcher_mean_spd Works", {
  X <- array(0, dim = c(3,3,5))
  set.seed(12345)
  for(i in 1:5) {
    X[,,i] <- MGLMRiem::randspd(n = 3)
  }  

  expect_equal(class(karcher_mean_spd(X, niter=5)), "matrix")
})
