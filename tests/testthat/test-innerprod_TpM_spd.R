test_that("innerprod_TpM_spd", {
  expect_equal(innerprod_TpM_spd(U = diag(rep(1,3)),
                                 V = diag(rep(1,3)),
                                 P = diag(rep(1,3))), 3)
  
  P= randspd(3)
  U= randspd(3)
  V= randspd(3)
  expect_equal(length(innerprod_TpM_spd(U = U,V = V,P = P)), 1)
})
