test_that("results are deterministic for the same input", {
  t <- seq(0, 2, length.out = 1200)
  x <- cos(2*pi*3*t) + cos(2*pi*11*t)

  r1 <- rfif(x)
  r2 <- rfif(x)

  expect_equal(r1$nimf, r2$nimf)
  expect_equal(r1$imfs, r2$imfs, tolerance = 1e-10)
  expect_equal(r1$residual, r2$residual, tolerance = 1e-10)
})
