\name{rfif}
\alias{rfif}
\title{Fast Iterative Filtering (C backend)}
\description{
Decompose a univariate signal into Intrinsic Mode Functions (IMFs) using the vendored C implementation
of Fast Iterative Filtering (FIF).
}
\usage{
rfif(x)
}
\arguments{
\item{x}{Numeric vector (signal).}
}
\value{
A list with components:
\itemize{
  \item \code{imfs}: numeric matrix with one IMF per row (rows = IMFs; columns = time).
  \item \code{residual}: numeric vector of length \code{length(x)}.
  \item \code{nimf}: integer number of IMFs.
}
}
\examples{
t <- seq(0, 1, length.out = 1024)
x <- sin(2*pi*5*t) + 0.5*sin(2*pi*20*t)
res <- rfif(x)
recon <- if (res$nimf > 0) colSums(res$imfs) + res$residual else res$residual
max(abs(x - recon))
}
