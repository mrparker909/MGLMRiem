test_that("isspd Works", {
  # testing for SPD T/F
  expect_equal(isspd(diag(rep(1, 3))), 1)
  expect_equal(isspd(diag(rep(0, 3))), 0)
  
  # testing for symmetric T/F
  expect_equal(issym(diag(rep(1, 3))), T)
  expect_equal(issym(matrix(c(1,2,2,2), nrow=2)), T)
})
