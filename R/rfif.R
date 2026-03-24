#' Fast Iterative Filtering (C backend using R's internal FFT)
#'
#' @param x Numeric vector (signal).
#' @return A list with components:
#' \describe{
#'   \item{imfs}{Numeric matrix (nimf x length(x)) of intrinsic mode functions (IMFs).}
#'   \item{residual}{Numeric vector, the residual component.}
#'   \item{nimf}{Integer number of IMFs returned.}
#' }
#' @export
rfif <- function(x) {
  x <- as.numeric(x)

  # Pad to the next power of two for robustness across toolchains, then truncate outputs.
  .rfif_nextpow2 <- function(n) {
    if (n <= 1L) return(1L)
    as.integer(2^(ceiling(log2(as.double(n)))))
  }

  n0 <- length(x)
  N  <- .rfif_nextpow2(n0)
  xpad <- if (N != n0) c(x, rep(0, N - n0)) else x

  res <- tryCatch(
    .Call("r_fifc_run", xpad, PACKAGE = "RFIF"),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("no IMFs", msg, ignore.case = TRUE)) {
        list(imfs = matrix(0, nrow = 0, ncol = length(xpad)), residual = xpad, nimf = 0L)
      } else {
        stop(e)
      }
    }
  )

# Truncate (or guard against) any padding back to the original length
  if (!is.null(res$imfs) && is.matrix(res$imfs) && ncol(res$imfs) != n0) {
    res$imfs <- res$imfs[, seq_len(n0), drop = FALSE]
  }
  if (!is.null(res$residual) && length(res$residual) != n0) {
    res$residual <- res$residual[seq_len(n0)]
  }

  res
}
