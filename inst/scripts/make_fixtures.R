# Generates the reference fixtures in tests/testthat/fixtures/.
#
# Run from the package root after installing the package:
#     Rscript inst/scripts/make_fixtures.R
#
# The script writes three input cases as headerless CSV files with 17
# significant digits, so that R, Python and MATLAB read identical double
# values:
#   xy.csv          60 x 80 cross-omics matrix with two planted blocks
#   x.csv           50 x 50 symmetric matrix with two planted communities
#   three_step_P.csv, three_step_R.csv
#                   70 x 70 association matrices (m = 30, n = 40) computed
#                   by triger_association() from simulated data in which
#                   X features 1-8 load on factor f1, X features 9-16 load on
#                   factor f2, and Y features 1-20 load on f1 + f2
#
# It then writes the expected outputs of the TriGer functions on these
# inputs (see write_expected() below). inst/scripts/make_fixtures.m writes
# the same expected-output files from the MATLAB implementation; the files in
# tests/testthat/fixtures/ are those of the last script run, as recorded in
# tests/testthat/fixtures/README.md.
#
# The Python package (TriGer-PyPI) uses copies of these files in
# tests/fixtures/.

library(TriGer)

out <- file.path("tests", "testthat", "fixtures")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

write_matrix <- function(M, file) {
    lines <- apply(M, 1L, function(r) {
        paste(sprintf("%.17g", r), collapse = ",")
    })
    writeLines(lines, file.path(out, file))
}

set.seed(20261121)

XY <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
XY[6:20, 11:30] <- XY[6:20, 11:30] + 4
XY[31:40, 51:65] <- XY[31:40, 51:65] + 3
write_matrix(XY, "xy.csv")

S <- matrix(runif(50 * 50, -0.2, 0.2), 50, 50)
S[1:12, 1:12] <- S[1:12, 1:12] + 0.9
S[25:32, 25:32] <- S[25:32, 25:32] + 0.7
S <- (S + t(S)) / 2
write_matrix(S, "x.csv")

N <- 100
m <- 30
n <- 40
f1 <- rnorm(N)
f2 <- rnorm(N)
X <- matrix(rnorm(N * m), N)
Y <- matrix(rnorm(N * n), N)
X[, 1:8] <- X[, 1:8] + 1.2 * f1
X[, 9:16] <- X[, 9:16] + 1.2 * f2
Y[, 1:20] <- Y[, 1:20] + 0.9 * (f1 + f2)
assoc <- triger_association(X, Y)
write_matrix(unname(assoc$P), "three_step_P.csv")
write_matrix(unname(assoc$R), "three_step_R.csv")

# Expected outputs, with the settings used by the tests:
#   xy_blocks.csv, xy_density.csv
#       greedy_peeling_XY_all(xy, lam = 1.5, max_blocks = 4)
#   x_blocks.csv, x_density.csv
#       greedy_peeling_X_all(x, lam = 1.5)
#   three_step_blocks.csv
#       greedy_3step(P, R, m = 30, n = 40, lambda1 = 1.1, lambda2 = 1.5,
#           lambda3 = 1.5, X_size = c(5, 10), Y_size = c(5, 15),
#           max_blocks = 4)
# Block files have columns block, side and index (1-based), with the indices
# of each block side in output order. Density files have columns block and
# density.
write_blocks <- function(blocks, sides, file) {
    rows <- lapply(seq_along(blocks), function(k) {
        do.call(rbind, lapply(sides, function(s) {
            data.frame(block = k, side = s, index = blocks[[k]][[s]])
        }))
    })
    utils::write.csv(do.call(rbind, rows), file.path(out, file),
        row.names = FALSE, quote = FALSE)
}

write_density <- function(blocks, file) {
    d <- data.frame(block = seq_along(blocks),
        density = sprintf("%.17g",
            vapply(blocks, `[[`, numeric(1), "density")))
    utils::write.csv(d, file.path(out, file), row.names = FALSE, quote = FALSE)
}

xy_blocks <- greedy_peeling_XY_all(XY, lam = 1.5, max_blocks = 4)
write_blocks(xy_blocks, c("X", "Y"), "xy_blocks.csv")
write_density(xy_blocks, "xy_density.csv")

x_comms <- greedy_peeling_X_all(S, lam = 1.5)
write_blocks(x_comms, "X", "x_blocks.csv")
write_density(x_comms, "x_density.csv")

modules <- greedy_3step(assoc$P, assoc$R, m, n, lambda1 = 1.1,
    lambda2 = 1.5, lambda3 = 1.5, X_size = c(5, 10), Y_size = c(5, 15),
    max_blocks = 4)
write_blocks(modules, c("X", "Y"), "three_step_blocks.csv")
