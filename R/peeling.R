#' Peel one dense block from a cross-omics matrix
#'
#' Greedily removes the X feature (row) or Y feature (column) with the
#' weakest remaining association and keeps the intermediate submatrix with
#' the highest penalized density, `sum / (|X| * |Y|)^(lam / 2)`. A row is
#' removed when its sum, scaled by the current ratio `|X| / |Y|`, is at most
#' the smallest column sum; otherwise a column is removed.
#'
#' @param Wp Non-negative matrix of cross-omics associations (X features in
#'   rows, Y features in columns), for example -log10 p-values.
#' @param lam Size penalty, typically in \[1, 2\]. Higher values give smaller,
#'   denser blocks.
#' @return A list with `density` (the penalized density of the block), `X`
#'   (row indices) and `Y` (column indices).
#' @examples
#' set.seed(1)
#' W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
#' W[6:20, 11:30] <- W[6:20, 11:30] + 4
#' blk <- greedy_peeling_XY_one(W, 1.5)
#' range(blk$X)
#' range(blk$Y)
#' @export
greedy_peeling_XY_one <- function(Wp, lam) {
    Wp <- .as_numeric_matrix(Wp, "Wp")
    .check_number(lam, "lam")
    m <- nrow(Wp)
    n <- ncol(Wp)
    C <- colSums(Wp)
    R <- rowSums(Wp)
    rem_X <- logical(m)
    rem_Y <- logical(n)
    best_rem_X <- rem_X
    best_rem_Y <- rem_Y
    len_X <- m
    len_Y <- n
    W_density <- 0
    for (it in seq_len(m + n)) {
        Cc <- C
        Cc[rem_Y] <- Inf
        Rr <- R
        Rr[rem_X] <- Inf
        IndT <- which.min(Cc)
        dT <- Cc[IndT]
        IndS <- which.min(Rr)
        dS <- Rr[IndS]
        if (!is.finite(dT) || !is.finite(dS)) break
        if ((len_X / len_Y) * dS <= dT) {
            rem_X[IndS] <- TRUE
            C <- C - Wp[IndS, ]
            C[rem_Y] <- 0
            Wp_sum <- sum(C)
            len_X <- len_X - 1L
        } else {
            rem_Y[IndT] <- TRUE
            R <- R - Wp[, IndT]
            R[rem_X] <- 0
            Wp_sum <- sum(R)
            len_Y <- len_Y - 1L
        }
        if (len_X <= 0L || len_Y <= 0L) break
        score <- Wp_sum / (len_X * len_Y)^(lam / 2)
        if (is.finite(score) && score > W_density) {
            W_density <- score
            best_rem_X <- rem_X
            best_rem_Y <- rem_Y
        }
    }
    list(density = W_density, X = which(!best_rem_X), Y = which(!best_rem_Y))
}


#' Peel dense blocks from a cross-omics matrix repeatedly
#'
#' Applies [greedy_peeling_XY_one()] repeatedly, removing the rows and
#' columns of each accepted block before the next round. A block is accepted
#' while its penalized density exceeds that of the full matrix, up to
#' `max_blocks` blocks.
#'
#' @param XY Non-negative cross-omics association matrix.
#' @inheritParams greedy_peeling_XY_one
#' @param max_blocks Maximum number of blocks returned.
#' @return A list of blocks, each `list(X, Y, density)`, with indices into
#'   the rows (`X`) and columns (`Y`) of `XY`.
#' @examples
#' set.seed(1)
#' W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
#' W[6:20, 11:30] <- W[6:20, 11:30] + 4
#' blocks <- greedy_peeling_XY_all(W, 1.5)
#' lengths(blocks[[1]][c("X", "Y")])
#' @export
greedy_peeling_XY_all <- function(XY, lam = 1.5, max_blocks = 5) {
    XY <- .as_numeric_matrix(XY, "XY")
    .check_number(lam, "lam")
    .check_count(max_blocks, "max_blocks")
    org_density <- sum(XY) / (nrow(XY) * ncol(XY))^(lam / 2)
    remain_X <- seq_len(nrow(XY))
    remain_Y <- seq_len(ncol(XY))
    result <- list()
    while (length(result) < max_blocks &&
        length(remain_X) > 0L && length(remain_Y) > 0L) {
        pk <- greedy_peeling_XY_one(XY[remain_X, remain_Y, drop = FALSE], lam)
        if (pk$density <= org_density ||
            length(pk$X) == 0L || length(pk$Y) == 0L) {
            break
        }
        blk <- list(X = remain_X[pk$X], Y = remain_Y[pk$Y],
            density = pk$density)
        result[[length(result) + 1L]] <- blk
        remain_X <- setdiff(remain_X, blk$X)
        remain_Y <- setdiff(remain_Y, blk$Y)
    }
    result
}


#' Peel one dense community from a within-omics matrix
#'
#' Greedily removes the feature with the smallest remaining sum and keeps the
#' intermediate subset with the highest penalized density, `sum / size^lam`.
#'
#' @param Wp Symmetric matrix of within-omics associations, for example
#'   signed correlations.
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
    Wp <- .as_numeric_matrix(Wp, "Wp")
    .check_symmetric(Wp, "Wp")
    .check_number(lam, "lam")
    m <- nrow(Wp)
    C <- colSums(Wp)
    rem <- logical(m)
    best_rem <- rem
    W_density <- 0
    for (i in seq_len(m)) {
        Cc <- C
        Cc[rem] <- Inf
        Ind <- which.min(Cc)
        rem[Ind] <- TRUE
        C <- C - Wp[Ind, ]
        C[rem] <- 0
        if (i == m) break
        score <- sum(C) / (m - i)^lam
        if (is.finite(score) && score > W_density) {
            W_density <- score
            best_rem <- rem
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
#' @param min_size Smallest community returned. Communities with fewer
#'   features are peeled but not returned. The default, 1, returns every
#'   community.
#' @return A list of communities, each `list(X, density)`, with indices into
#'   the rows of `Wp`.
#' @examples
#' set.seed(1)
#' S <- matrix(runif(50 * 50, -0.2, 0.2), 50, 50)
#' S[1:12, 1:12] <- S[1:12, 1:12] + 0.9
#' S <- (S + t(S)) / 2
#' comms <- greedy_peeling_X_all(S, 1.5)
#' comms[[1]]$X
#' @export
greedy_peeling_X_all <- function(Wp, lam, min_size = 1) {
    Wp <- .as_numeric_matrix(Wp, "Wp")
    .check_symmetric(Wp, "Wp")
    .check_number(lam, "lam")
    .check_count(min_size, "min_size")
    org_density <- sum(Wp) / nrow(Wp)^lam
    remain <- seq_len(nrow(Wp))
    result <- list()
    while (length(remain) > 0L) {
        pk <- greedy_peeling_X_one(Wp[remain, remain, drop = FALSE], lam)
        if (pk$density <= org_density || length(pk$X) == 0L) break
        comm <- list(X = remain[pk$X], density = pk$density)
        result[[length(result) + 1L]] <- comm
        remain <- setdiff(remain, comm$X)
    }
    Filter(function(cm) length(cm$X) >= min_size, result)
}


#' Order features within detected blocks
#'
#' `sort_result()` orders the features of each block by their within-block
#' sums in `Wp`; for a cross-omics `Wp`, rows and columns are ordered by the
#' absolute values of their sums. `sort_result_XandY()` orders X features by
#' their sums within `X[blk$X, blk$X]` and Y features by their sums within
#' `Y[blk$Y, blk$Y]`.
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
#' head(blocks[[1]]$X)
#' @export
sort_result <- function(Wp, blocks, symmetric = NULL) {
    Wp <- as.matrix(Wp)
    if (is.null(symmetric)) symmetric <- nrow(Wp) == ncol(Wp)
    lapply(blocks, function(blk) {
        if (symmetric) {
            idx <- blk$X
            sums <- rowSums(Wp[idx, idx, drop = FALSE])
            blk$X <- idx[order(sums, decreasing = TRUE)]
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
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    lapply(blocks, function(blk) {
        sx <- colSums(X[blk$X, blk$X, drop = FALSE])
        sy <- colSums(Y[blk$Y, blk$Y, drop = FALSE])
        blk$X <- blk$X[order(sx, decreasing = TRUE)]
        blk$Y <- blk$Y[order(sy, decreasing = TRUE)]
        blk
    })
}


#' Detect cross-omics modules from association matrices
#'
#' The TriGer pipeline on precomputed matrices. Step 1 peels dense blocks
#' from the cross-omics part of `P` with [greedy_peeling_XY_all()]. Blocks
#' with `X_size[1]` or fewer X features, or `Y_size[1]` or fewer Y features,
#' are dropped. In steps 2 and 3, the X and Y sides of blocks with more than
#' `X_size[2]` or `Y_size[2]` features are refined by within-omics peeling of
#' the signed correlations in `R`: the side is replaced by the communities
#' with at least `min_community_size` features, each ordered by within-
#' community sums. [triger()] computes `P` and `R` from data and calls this
#' function.
#'
#' @param P `(m + n) x (m + n)` non-negative association matrix over the X
#'   and Y features (for example -log10 p-values of their correlations); only
#'   `P[1:m, (m + 1):(m + n)]` is used.
#' @param R `(m + n) x (m + n)` signed correlation matrix over the X and Y
#'   features, used for refinement.
#' @param m,n Numbers of X and Y features.
#' @param lambda1 Cross-omics size penalty (step 1).
#' @param lambda2,lambda3 Within-omics size penalties for X and Y (steps 2
#'   and 3).
#' @param X_size,Y_size `c(min, max)`: blocks with `min` or fewer features on
#'   either side are dropped, and sides with more than `max` features are
#'   refined.
#' @param min_community_size Smallest within-omics community kept during
#'   refinement.
#' @param keep_unrefined_if_empty If refinement keeps no community, keep that
#'   side of the block unrefined (`TRUE`) or empty it (`FALSE`).
#' @param max_blocks Maximum number of blocks found in step 1.
#' @return A list of modules, each `list(X, Y, density)`, with `X` in `1:m`,
#'   `Y` in `1:n`, and `density` the penalized density of the step-1 block
#'   before refinement.
#' @examples
#' set.seed(1)
#' N <- 150
#' m <- 80
#' n <- 120
#' f <- rnorm(N)
#' X <- matrix(rnorm(N * m), N)
#' Y <- matrix(rnorm(N * n), N)
#' X[, 1:25] <- X[, 1:25] + 1.5 * f
#' Y[, 1:40] <- Y[, 1:40] + 1.5 * f
#' R <- cor(cbind(X, Y))
#' P <- abs(R)
#' diag(P) <- 0
#' modules <- greedy_3step(P, R, m, n, lambda1 = 1.1, lambda2 = 1.8,
#'     lambda3 = 1.8)
#' lengths(modules[[1]][c("X", "Y")])
#' @export
greedy_3step <- function(P, R, m, n, lambda1 = 1.5, lambda2 = 1.5,
    lambda3 = 1.5, X_size = c(10, 200), Y_size = c(10, 200),
    min_community_size = 2, keep_unrefined_if_empty = TRUE,
    max_blocks = 5) {
    P <- .as_numeric_matrix(P, "P")
    R <- .as_numeric_matrix(R, "R")
    .check_count(m, "m")
    .check_count(n, "n")
    if (nrow(P) != m + n || ncol(P) != m + n) {
        stop("`P` must be (m + n) x (m + n)", call. = FALSE)
    }
    if (nrow(R) != m + n || ncol(R) != m + n) {
        stop("`R` must be (m + n) x (m + n)", call. = FALSE)
    }
    .check_size(X_size, "X_size")
    .check_size(Y_size, "Y_size")
    .check_count(min_community_size, "min_community_size")

    XY <- P[seq_len(m), m + seq_len(n), drop = FALSE]
    blocks <- greedy_peeling_XY_all(XY, lambda1, max_blocks)
    blocks <- Filter(function(b) {
        length(b$X) > X_size[1] && length(b$Y) > Y_size[1]
    }, blocks)

    refine <- function(idx, offset, lam) {
        sub <- R[offset + idx, offset + idx, drop = FALSE]
        sub <- (sub + t(sub)) / 2
        comm <- greedy_peeling_X_all(sub, lam, min_size = min_community_size)
        if (length(comm) == 0L) {
            return(if (keep_unrefined_if_empty) idx else integer(0))
        }
        comm <- sort_result(sub, comm, symmetric = TRUE)
        idx[unlist(lapply(comm, `[[`, "X"))]
    }

    lapply(blocks, function(b) {
        if (length(b$X) > X_size[2]) b$X <- refine(b$X, 0L, lambda2)
        if (length(b$Y) > Y_size[2]) b$Y <- refine(b$Y, m, lambda3)
        b
    })
}
