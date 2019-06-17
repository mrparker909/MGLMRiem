test_that("paralleltranslateAtoB_spd Works", {
  a <- repmat(diag(rep(1,3)),2)
  b <- repmat(diag(rep(2,3)),2)
  w <- repmat(diag(rep(1,3)),2)
  
  paralleltranslateAtoB_spd(a,b,w)
  
  expect_equal(dist_M_spd(X = diag(rep(1,3)), Y = diag(rep(1,3))), 0)
  
  
})
