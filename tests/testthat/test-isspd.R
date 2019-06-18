test_that("isspd Works", {
  # testing for SPD T/F
  expect_equal(isspd(diag(rep(1, 3))), T)
  expect_equal(isspd(diag(rep(0, 3))), F)
  
  # testing for symmetric T/F
  expect_equal(issym(diag(rep(1, 3))), T)
  expect_equal(issym(matrix(c(1,2,2,2), nrow=2)), T)
  
  M <- matrix(c(31.51896, 26.08810, 24.34587,
                26.08810, 22.02744, 19.86721,
                24.34587, 19.86721, 19.79121),
              byrow = T, nrow=3)
  expect_equal(issym(M), T)
  expect_equal(isspd(M), T)
  
})
