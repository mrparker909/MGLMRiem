test_that("frechet_variance_spd", {
  set.seed(123)
  Y = randspd_FAST(3, NUM = 10)
  YK=karcher_mean_spd(Y,niter=100)
  a=frechet_variance_spd(Y, maxiter = 100)
  b=frechet_variance_spd(Y,YK)
  expect_equal(a,b)
})
