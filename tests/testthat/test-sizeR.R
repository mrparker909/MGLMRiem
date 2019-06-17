test_that("sizeR Works", {
  X <- array(0, dim = c(3,3,5))
  
  expect_equal(sizeR(M = X), c(3,3,5))
})
