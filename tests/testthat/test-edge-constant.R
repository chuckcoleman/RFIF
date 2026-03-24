test_that("constant signal does not crash and reconstructs", {
  x <- rep(5, 800)
  res <- rfif(x)

  recon <- colSums(res$imfs) + res$residual
  expect_lt(max(abs(recon - x)), 1e-8)
})
