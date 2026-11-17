#' Peel one dense block from a cross-omics matrix
#'
#' Greedily removes the X feature (row) or Y feature (column) with the
#' weakest remaining association and keeps the intermediate submatrix with the
#' highest penalized density, `sum / (|X| * |Y|)^(lam / 2)`.
#'
#' @param Wp Non-negative matrix of cross-omics associations (X features in
#'   rows, Y features in columns), for example -log10 p-values.
#' @param lam Size penalty in \[1, 2\]. Higher values give smaller, denser blocks.
#' @return A list with `density`, `X` (row indices) and `Y` (column indices).
#' @examples
#' set.seed(1)
#' W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
#' W[6:20, 11:30] <- W[6:20, 11:30] + 4
#' blk <- greedy_peeling_XY_one(W, 1.5)
#' range(blk$X); range(blk$Y)
#' @export
greedy_peeling_XY_one <- function(Wp, lam) {
  Wp <- as.matrix(Wp); storage.mode(Wp) <- "double"
  m <- nrow(Wp); n <- ncol(Wp)
  C <- colSums(Wp); R <- rowSums(Wp)
  rem_X <- logical(m); rem_Y <- logical(n)
  best_rem_X <- logical(m); best_rem_Y <- logical(n)
  len_X <- m; len_Y <- n
  W_density <- 0
  for (it in seq_len(m + n)) {
    Cc <- C; Cc[rem_Y] <- Inf
    Rr <- R; Rr[rem_X] <- Inf
    IndT <- which.min(Cc); dT <- Cc[IndT]
    IndS <- which.min(Rr); dS <- Rr[IndS]
    if (!is.finite(dT) || !is.finite(dS)) break
    cratio <- if (len_Y > 0) len_X / len_Y else Inf
    if (cratio * dS <= dT) {
      rem_X[IndS] <- TRUE
      C <- C - Wp[IndS, ]; C[rem_Y] <- 0
      Wp_sum <- sum(C); len_X <- len_X - 1L
    } else {
      rem_Y[IndT] <- TRUE
      R <- R - Wp[, IndT]; R[rem_X] <- 0
      Wp_sum <- sum(R); len_Y <- len_Y - 1L
    }
    if (len_X <= 0 || len_Y <= 0) break
    score <- Wp_sum / (len_X * len_Y)^(lam / 2)
    if (is.finite(score) && score > W_density) {
      W_density <- score; best_rem_X <- rem_X; best_rem_Y <- rem_Y
    }
  }
  list(density = W_density, X = which(!best_rem_X), Y = which(!best_rem_Y))
}


#' Peel dense blocks from a cross-omics matrix repeatedly
#'
#' Applies [greedy_peeling_XY_one()] repeatedly, removing each accepted block
#' before the next round. A block is accepted while its penalized density
#' exceeds that of the full matrix.
#'
#' @param XY Non-negative cross-omics association matrix.
#' @inheritParams greedy_peeling_XY_one
#' @param max_round Number of rounds; at most `max_round - 1` blocks are returned.
#' @return A list of blocks, each `list(X, Y, density)` with indices into `XY`.
#' @examples
#' set.seed(1)
#' W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
#' W[6:20, 11:30] <- W[6:20, 11:30] + 4
#' blocks <- greedy_peeling_XY_all(W, 1.5)
#' lengths(blocks[[1]][c("X", "Y")])
#' @export
greedy_peeling_XY_all <- function(XY, lam = 1.5, max_round = 5) {
  XY <- as.matrix(XY); storage.mode(XY) <- "double"
  m <- nrow(XY); n <- ncol(XY)
  org_density <- sum(XY) / (m * n)^(lam / 2)
  remain_X <- seq_len(m); remain_Y <- seq_len(n)
  Wp <- XY; result <- list()
  for (i in seq_len(max_round)) {
    if (length(remain_X) < 1L || length(remain_Y) < 1L) break
    pk <- greedy_peeling_XY_one(Wp, lam)
    if (pk$density > org_density && i < max_round &&
        length(pk$X) > 0 && length(pk$Y) > 0) {
      result[[length(result) + 1L]] <- list(X = sort(remain_X[pk$X]),
                                            Y = sort(remain_Y[pk$Y]),
                                            density = pk$density)
      remain_X <- setdiff(seq_len(m), unlist(lapply(result, `[[`, "X")))
      remain_Y <- setdiff(seq_len(n), unlist(lapply(result, `[[`, "Y")))
      if (length(remain_X) == 0L || length(remain_Y) == 0L) break
      Wp <- XY[remain_X, remain_Y, drop = FALSE]
    } else break
  }
  result
}


#' Peel one dense community from a within-omics matrix
#'
#' Greedily removes the feature with the smallest remaining sum and keeps the
#' intermediate subset with the highest penalized density, `sum / size^lam`.
#'
#' @param Wp Symmetric matrix of within-omics associations.
#' @inheritParams greedy_peeling_XY_one
#' @return A list with `density` and `X` (indices of the community).
#' @examples
#' set.seed(1)
#' S <- matrix(runif(50 * 50, -0.2, 0.2), 50, 50)
#' S[1:12, 1:12] <- S[1:12, 1:12] + 0.9
#' S <- (S + t(S)) / 2
#' greedy_peeling_X_one(S, 1.5)$X
#' @export
greedy_peeling_X_one <- function(Wp, lam) {
  Wp <- as.matrix(Wp); storage.mode(Wp) <- "double"
  m <- nrow(Wp)
  C <- colSums(Wp)
  rem <- logical(m); best_rem <- logical(m)
  W_density <- 0
  for (i in seq_len(m)) {
    Cc <- C; Cc[rem] <- Inf
    Ind <- which.min(Cc)
    rem[Ind] <- TRUE
    C <- C - Wp[Ind, ]; C[rem] <- 0
    Wp_sum <- sum(C)
    denom <- (m - i)^lam
    score <- if (denom > 0) Wp_sum / denom else Inf
    if (is.finite(score) && score > W_density) {
      W_density <- score; best_rem <- rem
    }
  }
  list(density = W_density, X = which(!best_rem))
}


#' Peel dense communities from a within-omics matrix repeatedly
#'
#' Applies [greedy_peeling_X_one()] repeatedly, removing each accepted
#' community before the next round, while the community's penalized density
#' exceeds that of the full matrix.
#'
#' @inheritParams greedy_peeling_X_one
#' @return A list of communities, each `list(X, density)`.
#' @examples
#' set.seed(1)
#' S <- matrix(runif(50 * 50, -0.2, 0.2), 50, 50)
#' S[1:12, 1:12] <- S[1:12, 1:12] + 0.9
#' S <- (S + t(S)) / 2
#' comms <- greedy_peeling_X_all(S, 1.5)
#' comms[[1]]$X
#' @export
greedy_peeling_X_all <- function(Wp, lam) {
  Wp <- as.matrix(Wp); storage.mode(Wp) <- "double"
  if (!isTRUE(all.equal(Wp, t(Wp), tolerance = 1e-8, check.attributes = FALSE)))
    stop("`Wp` must be symmetric")
  m <- nrow(Wp)
  org_density <- sum(Wp) / m^lam
  remain <- seq_len(m); cur <- Wp; result <- list()
  repeat {
    pk <- greedy_peeling_X_one(cur, lam)
    if (pk$density > org_density && length(pk$X) > 0) {
      result[[length(result) + 1L]] <- list(X = remain[pk$X], density = pk$density)
      remain <- setdiff(seq_len(m), unlist(lapply(result, `[[`, "X")))
      if (length(remain) == 0L) break
      cur <- Wp[remain, remain, drop = FALSE]
    } else break
  }
  result
}


#' Order features within detected blocks
#'
#' `sort_result()` orders the features of each block by their within-block
#' sums in `Wp`. `sort_result_XandY()` orders X features by their sums within
#' `X[blk$X, blk$X]` and Y features by their sums within `Y[blk$Y, blk$Y]`.
#'
#' @param Wp The matrix the blocks were detected in.
#' @param blocks A list of blocks, as returned by [greedy_peeling_XY_all()],
#'   [greedy_peeling_X_all()] or [greedy_3step()].
#' @param symmetric Whether `Wp` is a within-omics matrix. Defaults to
#'   `TRUE` for square matrices.
#' @param X,Y Within-omics matrices for the X and Y features.
#' @return `blocks` with the features of each block reordered.
#' @examples
#' set.seed(1)
#' W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
#' W[6:20, 11:30] <- W[6:20, 11:30] + 4
#' blocks <- sort_result(W, greedy_peeling_XY_all(W, 1.5))
#' @export
sort_result <- function(Wp, blocks, symmetric = NULL) {
  Wp <- as.matrix(Wp)
  if (is.null(symmetric)) symmetric <- nrow(Wp) == ncol(Wp)
  lapply(blocks, function(blk) {
    if (symmetric) {
      idx <- blk$X
      blk$X <- idx[order(rowSums(Wp[idx, idx, drop = FALSE]), decreasing = TRUE)]
    } else {
      sub <- Wp[blk$X, blk$Y, drop = FALSE]
      blk$X <- blk$X[order(abs(rowSums(sub)), decreasing = TRUE)]
      blk$Y <- blk$Y[order(abs(colSums(sub)), decreasing = TRUE)]
    }
    blk
  })
}


#' @rdname sort_result
#' @export
sort_result_XandY <- function(X, Y, blocks) {
  X <- as.matrix(X); Y <- as.matrix(Y)
  lapply(blocks, function(blk) {
    blk$X <- blk$X[order(colSums(X[blk$X, blk$X, drop = FALSE]), decreasing = TRUE)]
    blk$Y <- blk$Y[order(colSums(Y[blk$Y, blk$Y, drop = FALSE]), decreasing = TRUE)]
    blk
  })
}


#' Detect cross-omics modules with TriGer
#'
#' The full TriGer pipeline. Step 1 peels dense blocks from the cross-omics
#' part of `P` with [greedy_peeling_XY_all()]. Blocks with too few X or Y
#' features are dropped. In steps 2 and 3, the X and Y sides of blocks larger
#' than `X_size[2]` or `Y_size[2]` are refined by within-omics peeling of the
#' signed correlations in `R`, keeping communities with at least
#' `min_community_size` features.
#'
#' @param P `(m + n) x (m + n)` non-negative association matrix over the X
#'   and Y features (for example -log10 p-values of their correlations); only
#'   `P[1:m, (m + 1):(m + n)]` is used.
#' @param R `(m + n) x (m + n)` signed correlation matrix over the X and Y
#'   features, used for refinement.
#' @param m,n Numbers of X and Y features.
#' @param lambda1 Cross-omics size penalty.
#' @param lambda2,lambda3 Within-omics size penalties for X and Y.
#' @param X_size,Y_size `c(min, max)`: blocks with `min` or fewer features on
#'   either side are dropped, and sides with more than `max` features are refined.
#' @param min_community_size Smallest within-omics community kept during refinement.
#' @param keep_unrefined_if_empty If refinement keeps no community, keep that
#'   side of the block unrefined (`TRUE`) or empty it (`FALSE`).
#' @param max_round Passed to [greedy_peeling_XY_all()].
#' @return A list of modules, each `list(X, Y, density)`, with `X` in `1:m`
#'   and `Y` in `1:n`.
#' @examples
#' set.seed(1)
#' N <- 150; m <- 80; n <- 120
#' f <- rnorm(N)
#' X <- matrix(rnorm(N * m), N); Y <- matrix(rnorm(N * n), N)
#' X[, 1:25] <- X[, 1:25] + 1.5 * f; Y[, 1:40] <- Y[, 1:40] + 1.5 * f
#' R <- cor(cbind(X, Y))
#' P <- abs(R); diag(P) <- 0
#' modules <- greedy_3step(P, R, m, n, lambda1 = 1.1, lambda2 = 1.8, lambda3 = 1.8)
#' lengths(modules[[1]][c("X", "Y")])
#' @export
greedy_3step <- function(P, R, m, n, lambda1 = 1.5, lambda2 = 1.5, lambda3 = 1.5,
                         X_size = c(10, 200), Y_size = c(10, 200),
                         min_community_size = 2, keep_unrefined_if_empty = TRUE,
                         max_round = 5) {
  P <- as.matrix(P); storage.mode(P) <- "double"
  R <- as.matrix(R); storage.mode(R) <- "double"
  if (nrow(P) != m + n || ncol(P) != m + n) stop("`P` must be (m + n) x (m + n)")
  if (nrow(R) != m + n || ncol(R) != m + n) stop("`R` must be (m + n) x (m + n)")

  blocks <- greedy_peeling_XY_all(P[seq_len(m), m + seq_len(n), drop = FALSE], lambda1, max_round)
  blocks <- Filter(function(b) length(b$X) > X_size[1] && length(b$Y) > Y_size[1], blocks)

  refine <- function(idx, sub, lam) {
    sub <- (sub + t(sub)) / 2
    comm <- greedy_peeling_X_all(sub, lam)
    comm <- Filter(function(cm) length(cm$X) >= min_community_size, comm)
    if (length(comm) == 0L) return(if (keep_unrefined_if_empty) idx else integer(0))
    comm <- sort_result(sub, comm, symmetric = TRUE)
    idx[unlist(lapply(comm, `[[`, "X"))]
  }

  lapply(blocks, function(b) {
    if (length(b$X) > X_size[2]) b$X <- refine(b$X, R[b$X, b$X, drop = FALSE], lambda2)
    if (length(b$Y) > Y_size[2]) b$Y <- refine(b$Y, R[m + b$Y, m + b$Y, drop = FALSE], lambda3)
    b
  })
}
