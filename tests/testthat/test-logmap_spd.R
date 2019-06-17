test_that("Log Map Works", {
  expect_equal(logmap_spd(P = diag(rep(1,3)), X = diag(rep(1,3))),diag(rep(0,3)))
  expect_equal(logmap_spd(diag(rep(exp(1),3)), expmap_spd(P = diag(rep(1,3)), X = diag(rep(1,3)))), diag(rep(0,3)))
})
