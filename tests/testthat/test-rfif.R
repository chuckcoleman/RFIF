test_that("rfif returns expected structure", {
  t <- seq(0, 1, length.out = 512)
  x <- sin(2*pi*5*t) + 0.5*sin(2*pi*20*t)

  res <- rfif(x)

  expect_type(res, "list")
  expect_true(all(c("imfs","residual","nimf") %in% names(res)))

  expect_type(res$nimf, "integer")
  expect_true(res$nimf >= 0)

  expect_true(is.matrix(res$imfs))
  expect_true(is.numeric(res$residual))

  expect_equal(ncol(res$imfs), length(x))
  expect_equal(length(res$residual), length(x))
  expect_equal(nrow(res$imfs), res$nimf)
})

test_that("rfif reconstructs the input (within tolerance)", {
  set.seed(1)
  n <- 512
  t <- seq(0, 1, length.out = n)
  x <- sin(2*pi*7*t) + 0.25*sin(2*pi*31*t) + 0.05*rnorm(n)

  res <- rfif(x)

  recon <- if (res$nimf > 0) colSums(res$imfs) + res$residual else res$residual
  expect_equal(as.numeric(recon), as.numeric(x), tolerance = 1e-5)
})

test_that("rfif is deterministic", {
  n <- 512
  t <- seq(0, 1, length.out = n)
  x <- sin(2*pi*9*t) + 0.1*sin(2*pi*40*t)

  r1 <- rfif(x)
  r2 <- rfif(x)

  expect_equal(r1$nimf, r2$nimf)
  expect_equal(r1$residual, r2$residual, tolerance = 0)
  expect_equal(r1$imfs, r2$imfs, tolerance = 0)
})

test_that("rfif handles constant signal", {
  x <- rep(3.14, 512)
  res <- rfif(x)
  recon <- if (res$nimf > 0) colSums(res$imfs) + res$residual else res$residual
  expect_equal(as.numeric(recon), as.numeric(x), tolerance = 1e-10)
  expect_true(all(is.finite(res$residual)))
  expect_true(all(is.finite(res$imfs)))
})
