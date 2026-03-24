# RFIF — Fast Iterative Filtering for R

RFIF provides an R interface to a C implementation of **Fast Iterative Filtering (FIF)** for decomposing a univariate signal into **intrinsic mode functions (IMFs)** plus a residual.

## Installation

Install from a source tarball:

```r
install.packages("RFIF_1.0.tar.gz", repos = NULL, type = "source")
```

## Fast FFT backend (FFTW3)

RFIF performs many FFTs during decomposition. To guarantee portability, the package ships with a **self-contained fallback FFT** that works everywhere, but it can be much slower.

For best performance, RFIF can use **FFTW3** when available.

### Backend selection

During installation:

- If FFTW3 is detected via `pkg-config` → RFIF enables the **fast FFTW backend**
- Otherwise → RFIF uses the **portable fallback FFT** (still fully functional)

When FFTW3 is found, installation prints:

```
Found fftw3 via pkg-config (enabling fast FFT).
```

## Installing FFTW3 (recommended)

### macOS (MacPorts)

```bash
sudo port install pkgconfig
sudo port install fftw-3
```

Verify:

```bash
pkg-config --modversion fftw3
```

### macOS (Homebrew)

```bash
brew install pkg-config fftw
```

### Linux (Ubuntu / Debian)

```bash
sudo apt-get install libfftw3-dev pkg-config
```

### Windows (Rtools + MSYS2)

1. Install Rtools: https://cran.r-project.org/bin/windows/Rtools/

2. Open the **Rtools MSYS2 shell**, then install FFTW and pkg-config:

```bash
pacman -S mingw-w64-x86_64-fftw
pacman -S mingw-w64-x86_64-pkg-config
```

Verify:

```bash
pkg-config --modversion fftw3
```

Then install RFIF from source in R as usual.

## Basic usage

```r
library(RFIF)

t <- seq(0, 1, length.out = 1000)
x <- sin(2*pi*5*t) + 0.5*sin(2*pi*20*t)

res <- rfif(x)

str(res)
```

Returned object:

- `imfs`: numeric matrix (rows = IMFs, columns = time)
- `residual`: numeric vector (same length as input)
- `nimf`: integer number of IMFs

## Reconstruction check

```r
recon <- if (res$nimf > 0) colSums(res$imfs) + res$residual else res$residual
max(abs(x - recon))
```

## Vignette

```r
browseVignettes("RFIF")
```
