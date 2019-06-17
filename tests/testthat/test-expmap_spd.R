test_that("Exponential Map Works", {
  expect_equal(expmap_spd(P = diag(rep(1,3)), X = diag(rep(1,3))), diag(rep(exp(1),3)))
  expect_equal(expmap_spd(P = diag(rep(1,3)), X = diag(rep(2,3))), diag(rep(exp(1)^2,3)))
})
