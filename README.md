# TriGer

[![R-CMD-check](https://github.com/MATRIX-Research-Lab/TriGer/actions/workflows/check.yml/badge.svg)](https://github.com/MATRIX-Research-Lab/TriGer/actions/workflows/check.yml)
[![bioRxiv](https://img.shields.io/badge/bioRxiv-10.64898%2F2026.09.10.749540-b31b1b)](https://doi.org/10.64898/2026.09.10.749540)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE.md)

| Platform | Repository | Installation |
|----------|------------|--------------|
| R | [MATRIX-Research-Lab/TriGer](https://github.com/MATRIX-Research-Lab/TriGer) (this repository) | `remotes::install_github("MATRIX-Research-Lab/TriGer")` |
| Python | [MATRIX-Research-Lab/TriGer-PyPI](https://github.com/MATRIX-Research-Lab/TriGer-PyPI) | `pip install triger` |
| MATLAB | [jy2941/TriGer](https://github.com/jy2941/TriGer) | Clone the repository and add `matlab/methods` and `matlab/visualization` to the MATLAB path |

The R and Python packages give identical results on the same input. The MATLAB
repository also contains the code for the simulations and the real-data
analyses of the paper.

R package for **TriGer**, from the paper *Structured cross-omics interaction
discovery with a triple-graph model*.

TriGer finds modules of features that are densely associated across two omics
layers (X and Y) and coherent within each layer. It peels dense blocks from the
cross-omics association matrix, then refines large blocks using the
within-omics correlation matrices.

## Installation

The package is being submitted to Bioconductor. Install the development
version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("MATRIX-Research-Lab/TriGer")
```

After acceptance, install it with:

```r
# install.packages("BiocManager")
BiocManager::install("TriGer")
```

## Quick start

```r
library(TriGer)

# x: samples x m features, y: samples x n features (matrices or data frames),
# or two SummarizedExperiment objects with features in rows
res <- triger(x, y, lambda1 = 1.1, lambda2 = 1.8, lambda3 = 1.8,
    X_size = c(10, 15), Y_size = c(10, 250))

# Module sizes and feature names
sapply(res$modules, function(b) c(X = length(b$X), Y = length(b$Y)))
res$modules[[1]]$X_names

# Cross- and within-omics correlations of the module features
m <- res$m
n <- res$n
plot3in1_v1(res$R[1:m, 1:m], res$R[m + 1:n, m + 1:n], res$R[1:m, m + 1:n],
    res$modules, line = TRUE)
```

`vignette("TriGer")` describes the method and works through a simulated
example.

## Functions

| Function | Description |
|----------|-------------|
| `triger()` | Full pipeline on two omics layers (matrices or `SummarizedExperiment`) |
| `triger_association()` | Signed correlation and -log10 p-value matrices over both layers |
| `greedy_3step()` | Full pipeline on precomputed association matrices |
| `greedy_peeling_XY_one()` / `greedy_peeling_XY_all()` | Dense blocks in a cross-omics (X-Y) matrix |
| `greedy_peeling_X_one()` / `greedy_peeling_X_all()` | Dense communities in a within-omics (symmetric) matrix |
| `sort_result()` / `sort_result_XandY()` | Order features within each block |
| `plot3in1_v1()` | X-Y heatmap with the X-X and Y-Y matrices as triangles |
| `plot_all_result()` | Heatmap of the detected blocks |

`lambda1` is the cross-omics size penalty and `lambda2` / `lambda3` are the
within-omics penalties for X and Y (λ₀, λ₁, λ₂ in the paper); higher values give
smaller, denser blocks. `max_blocks` is the maximum number of blocks found in
the cross-omics step (default 5). The CRC and IBD analyses of the paper, run
with the MATLAB code, correspond to `max_blocks = 4`.

## Citation

Yu J, Lin H, Chen S (2026). Structured cross-omics interaction discovery with a
triple-graph model. *bioRxiv*.
[doi:10.64898/2026.09.10.749540](https://doi.org/10.64898/2026.09.10.749540)

`citation("TriGer")` returns this reference in R.

## License

MIT
