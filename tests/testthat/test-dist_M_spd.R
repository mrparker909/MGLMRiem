test_that("dist_M_spd", {
  expect_equal(dist_M_spd(X = diag(rep(1,3)), Y = diag(rep(1,3))), 0)
  expect_equal(dist_M_spd(X = diag(rep(1,3)), Y = diag(rep(2,3))), 1.2, tolerance=10^-3)
})
