#' Fast Iterative Filtering (C backend using R's internal FFT)
#'
#' @param x Numeric vector (signal).
#' @return A list with components:
#' \describe{
#'   \item{imfs}{Numeric matrix (nimf x length(x)) of intrinsic mode functions.}
#'   \item{residual}{Numeric vector, x - colSums(imfs).}
#'   \item{nimf}{Integer number of IMFs returned.}
#' }
#' @export
rfif <- function(x) {
  x <- as.numeric(x)

  # Pad to next power-of-two for robustness on some toolchains/FFT backends,
  # then truncate outputs back to the original length.
  .rfif_nextpow2 <- function(n) {
    if (n <= 1L) return(1L)
    as.integer(2^(ceiling(log2(as.double(n)))))
  }
  n0 <- length(x)
  N <- .rfif_nextpow2(n0)
  if (N != n0) {
    x_pad <- c(x, rep(0, N - n0))
  } else {
    x_pad <- x
  }

  res <- .Call("r_fifc_run", x_pad, PACKAGE = "RFIF")

  if (N != n0) {
    if (!is.null(res$imfs) && is.matrix(res$imfs)) {
      res$imfs <- res$imfs[, seq_len(n0), drop = FALSE]
    }
    if (!is.null(res$residual)) {
      res$residual <- res$residual[seq_len(n0)]
    }
  }

  res
}
