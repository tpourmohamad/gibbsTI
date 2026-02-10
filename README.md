# gibbsTI <img src='man/figures/gibbsTI.png' align="right" height="240" />

### R Package: Calibrated Gibbs Posteriors for Tolerance Intervals

[![Project Status: Under Development](https://img.shields.io/badge/status-under--development-orange.svg)](https://www.repostatus.org/#active)

**gibbsTI** provides tools for constructing and calibrating Gibbs posteriors to create robust frequentist tolerance intervals. 

> ⚠️ **Note:** This package is currently under active development. The API is subject to change, and you may encounter bugs. Please report any issues or feedback via the [GitHub Issue Tracker](https://github.com/tpourmohamad/gibbsTI/issues).

---

## Installation

You can install the development version of **gibbsTI** from GitHub. 

### 1. Prerequisites
Since you are installing from source, ensure you have a functional build environment:
* **Windows:** [Rtools](https://cran.r-project.org/bin/windows/Rtools/)
* **macOS:** Xcode command line tools

### 2. Install via `remotes`
Run the following commands in your R console:

```R
# Install the remotes package if you don't have it
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}

# Install gibbsTI
remotes::install_github("tpourmohamad/gibbsTI")
