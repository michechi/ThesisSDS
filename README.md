# mfmsbm

Collapsed Gibbs Sampler for Mixture of Finite Mixtures of Stochastic Block
Models (MFM-SBM).

## Overview

**mfmsbm** implements Bayesian community detection in networks using the
MFM-SBM framework. It supports two priors on the number of components:

- **Zero-truncated Poisson** (Miller and Harrison, 2018)
- **Gnedin** (Gnedin, 2010)

The package provides functions for:

- Network generation from SBM (`generate_sbm()`)
- Collapsed Gibbs sampling (`cogibbs_poisson()`, `cogibbs_gnedin()`)
- V_n computation (`log_vn_miller()`, `log_vn_gnedin()`, `log_vn_approx()`)
- Post-processing: co-clustering matrices, edge probability estimation,
  posterior on K (`coclustering_matrix()`, `edge_probability()`, `posterior_k()`)
- Sensitivity analysis and convergence diagnostics

## Installation

```r
# Install from GitHub
# install.packages("devtools")
devtools::install_github("Mich9534/mfmsbm")
```

## Quick Start

```r
library(mfmsbm)

# Generate a network
set.seed(42)
sbm <- generate_sbm(n = 100, K = 3, p = 0.5, q = 0.1, type_network = "balanced")

# Run collapsed Gibbs sampler
logVn <- log_vn_miller(gamma_ = 1, n = 100, upto = 110, lambda = 1)
z_start <- c(sample(1:9, 9), sample(1:9, 91, replace = TRUE))
fit <- cogibbs_poisson(M = 500, K_start = 9, A = sbm$A0, n = 100,
                       a = 1, b = 1, gamma_ = 1, logVn_miller = logVn,
                       z_start = z_start)

# Post-processing
cc <- coclustering_matrix(fit$z_post[, 201:500])
```

## Citation

If you use this package, please cite:

> Chini, M. (2024). *Community Detection in Networks via Mixture of Finite
> Mixtures of Stochastic Block Models*. M.Sc. Thesis, University of Turin.

The methodology is based on:

> Geng, J., Bhattacharya, A., and Pati, D. (2019). "Probabilistic community
> detection with unknown number of communities."
> *Journal of the American Statistical Association*, 114(526), 893--905.
> doi:10.1080/01621459.2018.1458618

## License

MIT
