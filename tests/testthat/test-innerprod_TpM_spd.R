test_that("innerprod_TpM_spd", {
  expect_equal(innerprod_TpM_spd(U = diag(rep(1,3)),
                                 V = diag(rep(1,3)),
                                 P = diag(rep(1,3))), 3)
})
