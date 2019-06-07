test_that("multiplication works", {
  expect_equal(expmap_spd(P = diag(rep(1,3)), X = diag(rep(1,3))), diag(rep(exp(1),3)))
})
