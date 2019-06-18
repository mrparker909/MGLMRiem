test_that("Exponential Map Works", {
  expect_equal(expmap_spd(P = diag(rep(1,3)), X = diag(rep(1,3))), diag(rep(exp(1),3)))
  expect_equal(expmap_spd(P = diag(rep(1,3)), X = diag(rep(2,3))), diag(rep(exp(1)^2,3)))
  expect_equal(isspd(expmap_spd(P = diag(rep(1,3)), X = matrix(c(1,0,1,0,1,0,1,0,1), nrow=3))), T)
  
  expect_error(isspd(expmap_spd(P = diag(rep(1,3)), X = matrix(c(1,0,0,0,1,0,1,0,1), nrow=3))))
  
  X <- matrix(c(1,1,1,1,1,1,1,1,1), nrow=3)
  X_err <- matrix(c(-7.097457, -5.700037, -5.419913,
                    -5.700037, -4.967910, -4.351601,
                    -5.419913, -4.351601, -4.417031), nrow=3, byrow = T)
  M_err <- matrix(c(3.194412, 3.345391, 2.953779,
                    3.345391, 7.368579, 3.163411,
                    2.953779, 3.163411, 2.938806), nrow=3, byrow = T)
  M_goo <- matrix(c(4.626219, 3.615546, 3.918863,
                    3.615546, 4.502816, 2.622989,
                    3.918863, 2.622989, 4.622520), nrow=3, byrow = T)
  expect_equal(isspd(M_goo), T)
  expect_equal(isspd(M_err), T)
  expect_equal(isspd(X    ), F)
  expect_equal(isspd(X_err), F)
  expect_equal(isspd(expmap_spd(P = M_goo, X = X)), T)
  expect_equal(isspd(expmap_spd(P = M_err, X = X)), T)
  expect_equal(isspd(expmap_spd(P = M_goo, X = X_err)), T)
  expect_equal(isspd(expmap_spd(P = M_err, X = X_err)), T)
})

# # numerical error (P)
# [1,] 3.194412 3.345391 2.953779
# [2,] 3.345391 7.368579 3.163411
# [3,] 2.953779 3.163411 2.938806
# 
# # numerical error (X)
# [1,] -7.097457 -5.700037 -5.419913
# [2,] -5.700037 -4.967910 -4.351601
# [3,] -5.419913 -4.351601 -4.417031
#
# # no numerical error
# [1,] 4.626219 3.615546 3.918863
# [2,] 3.615546 4.502816 2.622989
# [3,] 3.918863 2.622989 4.622520