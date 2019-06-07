test_that("multiplication works", {
  expect_equal(logmap_spd(P = diag(rep(1,3)), X = diag(rep(1,3))),diag(rep(0,3)))
})
