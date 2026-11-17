planted <- function() {
  set.seed(1)
  W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
  W[6:20, 11:30] <- W[6:20, 11:30] + 4
  W
}

test_that("XY peeling recovers a planted block", {
  blk <- greedy_peeling_XY_one(planted(), 1.5)
  expect_true(all(6:20 %in% blk$X))
  expect_true(all(11:30 %in% blk$Y))
  blocks <- greedy_peeling_XY_all(planted(), 1.5)
  expect_equal(blocks[[1]]$X, 6:20)
  expect_equal(blocks[[1]]$Y, 11:30)
})

test_that("within-omics peeling finds a planted community", {
  set.seed(2)
  S <- matrix(runif(50 * 50, -0.2, 0.2), 50, 50)
  S[1:12, 1:12] <- S[1:12, 1:12] + 0.9
  S <- (S + t(S)) / 2
  expect_true(all(1:12 %in% greedy_peeling_X_all(S, 1.5)[[1]]$X))
  expect_error(greedy_peeling_X_all(matrix(1:4, 2), 1.5), "symmetric")
})

test_that("sorting keeps block members", {
  W <- planted()
  blocks <- greedy_peeling_XY_all(W, 1.5)
  s <- sort_result(W, blocks)
  expect_setequal(s[[1]]$X, blocks[[1]]$X)
  X <- cor(matrix(rnorm(30 * 60), 30)); Y <- cor(matrix(rnorm(30 * 80), 30))
  s2 <- sort_result_XandY(X, Y, blocks)
  expect_setequal(s2[[1]]$Y, blocks[[1]]$Y)
})

test_that("greedy_3step finds and refines the module", {
  set.seed(4)
  N <- 150; m <- 80; n <- 120
  f <- rnorm(N)
  X <- matrix(rnorm(N * m), N); Y <- matrix(rnorm(N * n), N)
  X[, 1:25] <- X[, 1:25] + 1.5 * f; Y[, 1:40] <- Y[, 1:40] + 1.5 * f
  R <- cor(cbind(X, Y)); P <- abs(R); diag(P) <- 0
  mods <- greedy_3step(P, R, m, n, 1.1, 1.8, 1.8, X_size = c(10, 20), Y_size = c(10, 30))
  expect_gte(length(mods), 1)
  expect_gte(length(intersect(mods[[1]]$X, 1:25)), 15)
  expect_true(all(mods[[1]]$Y <= n))
})
