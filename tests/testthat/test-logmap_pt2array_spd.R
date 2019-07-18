test_that("logmap_pt2array_spd Works", {
  X <- array(0, dim = c(3,3,5))
  set.seed(12345)
  for(i in 1:5) {
    X[,,i] <- MGLMRiem::randspd_FAST(n = 3)
  }  
  expect_equal(logmap_pt2array_spd(p=X[,,1],X)[,,1],diag(c(0,0,0)))
})
