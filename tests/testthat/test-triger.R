simulate_layers <- function(N = 100, m = 40, n = 30) {
    set.seed(5)
    f <- rnorm(N)
    x <- matrix(rnorm(N * m), N,
        dimnames = list(paste0("s", seq_len(N)), paste0("g", seq_len(m))))
    y <- matrix(rnorm(N * n), N,
        dimnames = list(paste0("s", seq_len(N)), paste0("m", seq_len(n))))
    x[, 1:12] <- x[, 1:12] + f
    y[, 1:8] <- y[, 1:8] + f
    list(x = x, y = y)
}

test_that("triger_association matches cor.test", {
    d <- simulate_layers()
    assoc <- triger_association(d$x, d$y)
    expect_equal(dim(assoc$R), c(70L, 70L))
    expect_equal(assoc$m, 40L)
    expect_equal(rownames(assoc$P)[c(1, 41)], c("g1", "m1"))
    expect_equal(unname(diag(assoc$P)), rep(0, 70))
    ct <- stats::cor.test(d$x[, 1], d$y[, 1])
    expect_equal(assoc$P["g1", "m1"], -log10(ct$p.value))
    expect_equal(assoc$R["g1", "m1"], unname(ct$estimate))
    sp <- triger_association(d$x, d$y, method = "spearman")
    expect_equal(sp$R["g2", "m3"],
        stats::cor(d$x[, 2], d$y[, 3], method = "spearman"))
})

test_that("triger_association rejects invalid input", {
    d <- simulate_layers()
    expect_error(triger_association(d$x, d$y[-1, ]), "same number of samples")
    x <- d$x
    x[, 3] <- 1
    expect_error(triger_association(x, d$y), "zero variance: g3")
    x[1, 1] <- NA
    expect_error(triger_association(x, d$y), "missing or infinite")
    expect_error(triger_association(d$x, d$y, method = "kendall"))
})

test_that("triger recovers the shared module with feature names", {
    d <- simulate_layers()
    res <- triger(d$x, d$y, lambda1 = 1.3, X_size = c(3, 50),
        Y_size = c(3, 50))
    expect_named(res, c("modules", "R", "P", "m", "n"))
    mod <- res$modules[[1]]
    expect_named(mod, c("X", "Y", "density", "X_names", "Y_names"))
    expect_gte(length(intersect(mod$X_names, paste0("g", 1:12))), 10)
    expect_gte(length(intersect(mod$Y_names, paste0("m", 1:8))), 6)
    expect_equal(mod$X_names, colnames(d$x)[mod$X])
})

test_that("triger matches samples by name and accepts SummarizedExperiment", {
    d <- simulate_layers()
    ref <- triger(d$x, d$y, lambda1 = 1.3, X_size = c(3, 50),
        Y_size = c(3, 50))
    shuffled <- triger(d$x, d$y[rev(rownames(d$y)), ], lambda1 = 1.3,
        X_size = c(3, 50), Y_size = c(3, 50))
    expect_equal(shuffled$modules, ref$modules)
    se_x <- SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = t(d$x)))
    se_y <- SummarizedExperiment::SummarizedExperiment(
        assays = list(abundance = t(d$y)))
    from_se <- triger(se_x, se_y, lambda1 = 1.3, X_size = c(3, 50),
        Y_size = c(3, 50), assay_x = "abundance")
    expect_equal(from_se$modules, ref$modules)
    y_bad <- d$y
    rownames(y_bad)[1] <- "other"
    expect_error(triger(d$x, y_bad), "sample names")
})

test_that("triger passes further arguments to greedy_3step", {
    d <- simulate_layers()
    res <- triger(d$x, d$y, lambda1 = 1.3, X_size = c(3, 50),
        Y_size = c(3, 50), max_blocks = 1)
    expect_length(res$modules, 1)
})
