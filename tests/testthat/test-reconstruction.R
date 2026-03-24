test_that("reconstruction identity holds: x = sum(IMFs) + residual", {
  set.seed(1)
  t <- seq(0, 1, length.out = 2000)
  x <- sin(2*pi*5*t) + 0.5*sin(2*pi*20*t)

  res <- rfif(x)

  expect_type(res, "list")
  expect_true(all(c("imfs", "residual", "nimf") %in% names(res)))
  expect_true(is.matrix(res$imfs))
  expect_equal(ncol(res$imfs), length(x))
  expect_equal(length(res$residual), length(x))

  recon <- colSums(res$imfs) + res$residual
  err <- max(abs(recon - x))

  # Tolerance chosen to be robust across platforms/BLAS/FFT backend details.
  expect_lt(err, 1e-6)
})
