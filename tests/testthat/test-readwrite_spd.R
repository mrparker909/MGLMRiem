test_that("read/write SPD works", {
  
  Y = randspd_FAST(3, NUM=10)
  writeSPD(Y)
  expect_equal(T, file.exists("mat.spd"))
  a = readSPD()
  expect_equal(Y, a$Y)
  file.remove("mat.spd")
})
