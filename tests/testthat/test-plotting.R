test_that("plots draw without error", {
    set.seed(1)
    W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
    W[6:20, 11:30] <- W[6:20, 11:30] + 4
    blocks <- sort_result(W, greedy_peeling_XY_all(W, 1.5))
    X <- cor(matrix(rnorm(40 * 60), 40))
    Y <- cor(matrix(rnorm(40 * 80), 40))
    f <- tempfile(fileext = ".pdf")
    grDevices::pdf(f, 5.6, 4.2)
    on.exit({
        grDevices::dev.off()
        unlink(f)
    })
    expect_silent(plot_all_result(W, blocks, full = TRUE))
    expect_silent(plot3in1_v1(X, Y, W, blocks, line = TRUE, step = c(10, 5),
        bar = TRUE))
    expect_length(matlab_jet(8), 8)
})
