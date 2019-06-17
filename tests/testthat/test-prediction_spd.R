test_that("prediction_spd Works", {
  p <- diag(rep(1,3)) # 3x3 SPD matrix
  V <- repmat(p, 5) # N=5 coefficients V
  X <- matrix(c(1,2,3,4,5,6,7,8,9,1,2,3,4,5,6)/10, nrow=3) # rows = Xi components (3), cols = number of covariates (5)
  
  expect_equal(class(prediction_spd(p = p,V = V, X = X)), "array")
})
