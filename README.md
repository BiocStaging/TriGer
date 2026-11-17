# TriGer

R package for **TriGer**, from the paper *Structured cross-omics interaction
discovery with a triple-graph model*.

TriGer finds modules of features that are densely associated across two omics
layers (X and Y) and coherent within each layer. It peels dense blocks from the
cross-omics association matrix, then refines large blocks using the
within-omics correlation matrices.

## Installation

```r
# install.packages("remotes")
remotes::install_github("jy2941/TriGer-R")
```

## Quick start

```r
library(TriGer)

# X: samples x m features, Y: samples x n features
N <- nrow(X); m <- ncol(X); n <- ncol(Y)
R <- cor(cbind(X, Y))                                  # signed correlations
t <- R * sqrt((N - 2) / pmax(1 - R^2, 1e-12))
P <- -log10(2 * pt(abs(t), N - 2, lower.tail = FALSE))  # -log10 p
diag(P) <- 0

modules <- greedy_3step(P, R, m, n, lambda1 = 1.1, lambda2 = 1.8, lambda3 = 1.8,
                        X_size = c(10, 15), Y_size = c(10, 250))
sapply(modules, function(b) c(X = length(b$X), Y = length(b$Y)))

plot3in1_v1(R[1:m, 1:m], R[m + 1:n, m + 1:n], R[1:m, m + 1:n], modules, line = TRUE)
```

## Functions

| Function | Description |
|----------|-------------|
| `greedy_3step()` | Full pipeline: cross-omics block detection, then within-omics refinement of large blocks |
| `greedy_peeling_XY_one()` / `greedy_peeling_XY_all()` | Dense blocks in a cross-omics (X-Y) matrix |
| `greedy_peeling_X_one()` / `greedy_peeling_X_all()` | Dense communities in a within-omics (symmetric) matrix |
| `sort_result()` / `sort_result_XandY()` | Order features within each block |
| `plot3in1_v1()` | X-Y heatmap with the X-X and Y-Y matrices as triangles |
| `plot_all_result()` | Heatmap of the detected blocks |

`lambda1` is the cross-omics size penalty and `lambda2` / `lambda3` are the
within-omics penalties for X and Y (λ₀, λ₁, λ₂ in the paper); higher values give
smaller, denser blocks.

## Related

- Python package: [TriGer-PyPI](https://github.com/jy2941/TriGer-PyPI)
- Code for the paper's simulations and analyses: [TriGer](https://github.com/jy2941/TriGer)

## License

MIT
