# Compares the package with the reference outputs described in
# tests/testthat/fixtures/README.md.

fixture <- function(file) test_path("fixtures", file)

read_matrix <- function(file) {
    unname(as.matrix(utils::read.csv(fixture(file), header = FALSE)))
}

read_blocks <- function(file, sides) {
    d <- utils::read.csv(fixture(file))
    lapply(split(d, d$block), function(b) {
        stats::setNames(lapply(sides, function(s) b$index[b$side == s]),
            sides)
    })
}

strip <- function(blocks, sides) {
    unname(lapply(blocks, function(b) b[sides]))
}

test_that("greedy_peeling_XY_all matches the reference output", {
    blocks <- greedy_peeling_XY_all(read_matrix("xy.csv"), 1.5,
        max_blocks = 4)
    expected <- read_blocks("xy_blocks.csv", c("X", "Y"))
    expect_equal(strip(blocks, c("X", "Y")), unname(expected))
    density <- utils::read.csv(fixture("xy_density.csv"))$density
    expect_equal(vapply(blocks, `[[`, numeric(1), "density"), density,
        tolerance = 1e-6)
})

test_that("greedy_peeling_X_all matches the reference output", {
    comms <- greedy_peeling_X_all(read_matrix("x.csv"), 1.5)
    expected <- read_blocks("x_blocks.csv", "X")
    expect_equal(strip(comms, "X"), unname(expected))
    density <- utils::read.csv(fixture("x_density.csv"))$density
    expect_equal(vapply(comms, `[[`, numeric(1), "density"), density,
        tolerance = 1e-6)
})

test_that("greedy_3step matches the reference output", {
    modules <- greedy_3step(read_matrix("three_step_P.csv"),
        read_matrix("three_step_R.csv"), m = 30, n = 40, lambda1 = 1.1,
        lambda2 = 1.5, lambda3 = 1.5, X_size = c(5, 10), Y_size = c(5, 15),
        max_blocks = 4)
    expected <- read_blocks("three_step_blocks.csv", c("X", "Y"))
    expect_equal(strip(modules, c("X", "Y")), unname(expected))
})
