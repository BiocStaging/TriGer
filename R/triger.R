#' Association matrices for two omics layers
#'
#' Computes the signed correlation matrix and the matrix of -log10 two-sided
#' p-values over the features of two omics layers measured on the same
#' samples. These are the inputs `R` and `P` of [greedy_3step()].
#'
#' The p-value of a correlation coefficient \eqn{r} estimated from \eqn{N}
#' samples is computed from \eqn{t = |r| \sqrt{(N - 2) / (1 - r^2)}} with
#' \eqn{N - 2} degrees of freedom, as in [stats::cor.test()] (for Spearman
#' correlations, the same t approximation). \eqn{r^2} is capped at
#' \eqn{1 - 10^{-12}} so that perfectly correlated features receive a large
#' finite value, and the diagonal of `P` is set to 0.
#'
#' @param x Numeric matrix of the first omics layer, samples in rows and
#'   features in columns.
#' @param y Numeric matrix of the second omics layer, with the same samples
#'   in the same order.
#' @param method Correlation coefficient: `"pearson"` or `"spearman"`.
#' @return A list with `R` and `P`, both `(m + n) x (m + n)` with the `m`
#'   features of `x` first and the `n` features of `y` second, and `m` and
#'   `n`. Row and column names are the column names of `x` and `y`, or
#'   `X1, X2, ...` and `Y1, Y2, ...` when these are missing.
#' @examples
#' set.seed(1)
#' x <- matrix(rnorm(50 * 5), 50)
#' y <- matrix(rnorm(50 * 4), 50)
#' assoc <- triger_association(x, y)
#' round(assoc$P[1:5, 6:9], 2)
#' @export
triger_association <- function(x, y, method = c("pearson", "spearman")) {
    method <- match.arg(method)
    x <- .as_numeric_matrix(x, "x")
    y <- .as_numeric_matrix(y, "y")
    if (nrow(x) != nrow(y)) {
        stop("`x` and `y` must have the same number of samples (rows)",
            call. = FALSE)
    }
    N <- nrow(x)
    if (N < 3L) stop("at least 3 samples are required", call. = FALSE)
    xy <- cbind(x, y)
    colnames(xy) <- c(.feature_names(x, "X"), .feature_names(y, "Y"))
    constant <- apply(xy, 2L, function(v) all(v == v[1L]))
    if (any(constant)) {
        stop("features with zero variance: ",
            paste(utils::head(colnames(xy)[constant], 5L), collapse = ", "),
            call. = FALSE)
    }
    R <- stats::cor(xy, method = method)
    t_stat <- abs(R) * sqrt((N - 2) / (1 - pmin(R^2, 1 - 1e-12)))
    log_p <- log(2) +
        stats::pt(t_stat, df = N - 2, lower.tail = FALSE, log.p = TRUE)
    P <- -log_p / log(10)
    diag(P) <- 0
    list(R = R, P = P, m = ncol(x), n = ncol(y))
}


#' Detect cross-omics modules with TriGer
#'
#' Runs the TriGer pipeline on two omics layers measured on the same
#' samples. [triger_association()] computes the cross- and within-omics
#' correlations and their -log10 p-values, and [greedy_3step()] finds blocks
#' of X and Y features that are densely associated across the two layers and
#' refines large blocks by within-omics peeling.
#'
#' Samples are matched by name when both layers carry sample names (row
#' names of a matrix, column names of a `SummarizedExperiment`); `y` is then
#' reordered to the sample order of `x`. Without sample names, the samples
#' must be in the same order.
#'
#' @param x,y The two omics layers: numeric matrices or data frames with
#'   samples in rows and features in columns, or
#'   [SummarizedExperiment::SummarizedExperiment] objects with features in
#'   rows and samples in columns.
#' @param lambda1 Cross-omics size penalty (step 1). Higher values give
#'   smaller, denser blocks.
#' @param lambda2,lambda3 Within-omics size penalties for X and Y (steps 2
#'   and 3).
#' @param X_size,Y_size `c(min, max)` passed to [greedy_3step()]: blocks with
#'   `min` or fewer features on either side are dropped, and sides with more
#'   than `max` features are refined.
#' @param method Correlation coefficient: `"pearson"` or `"spearman"`.
#' @param assay_x,assay_y Assay used when `x` or `y` is a
#'   `SummarizedExperiment`, by name or index.
#' @param ... Further arguments to [greedy_3step()]: `min_community_size`,
#'   `keep_unrefined_if_empty`, `max_blocks`.
#' @return A list with
#'   \describe{
#'     \item{`modules`}{A list of modules, each `list(X, Y, density, X_names,
#'       Y_names)`, where `X` and `Y` index the features of `x` and `y` and
#'       `density` is the penalized density of the step-1 block.}
#'     \item{`R`, `P`, `m`, `n`}{The output of [triger_association()].}
#'   }
#' @examples
#' set.seed(1)
#' N <- 100
#' f <- rnorm(N)
#' x <- matrix(rnorm(N * 40), N, dimnames = list(NULL, paste0("g", 1:40)))
#' y <- matrix(rnorm(N * 30), N, dimnames = list(NULL, paste0("m", 1:30)))
#' x[, 1:12] <- x[, 1:12] + f
#' y[, 1:8] <- y[, 1:8] + f
#' res <- triger(x, y, lambda1 = 1.3, X_size = c(3, 50), Y_size = c(3, 50))
#' res$modules[[1]]$X_names
#' res$modules[[1]]$Y_names
#' @export
triger <- function(x, y, lambda1 = 1.5, lambda2 = 1.5, lambda3 = 1.5,
    X_size = c(10, 200), Y_size = c(10, 200),
    method = c("pearson", "spearman"), assay_x = 1, assay_y = 1, ...) {
    x <- .samples_by_features(x, assay_x, "x")
    y <- .samples_by_features(y, assay_y, "y")
    y <- .match_samples(x, y)
    assoc <- triger_association(x, y, method = method)
    m <- assoc$m
    n <- assoc$n
    modules <- greedy_3step(assoc$P, assoc$R, m, n, lambda1 = lambda1,
        lambda2 = lambda2, lambda3 = lambda3, X_size = X_size,
        Y_size = Y_size, ...)
    x_names <- rownames(assoc$R)[seq_len(m)]
    y_names <- rownames(assoc$R)[m + seq_len(n)]
    modules <- lapply(modules, function(b) {
        b$X_names <- x_names[b$X]
        b$Y_names <- y_names[b$Y]
        b
    })
    c(list(modules = modules), assoc)
}


.samples_by_features <- function(obj, assay, arg) {
    if (methods::is(obj, "SummarizedExperiment")) {
        return(t(SummarizedExperiment::assay(obj, assay)))
    }
    .as_numeric_matrix(obj, arg)
}

.match_samples <- function(x, y) {
    sx <- rownames(x)
    sy <- rownames(y)
    if (is.null(sx) || is.null(sy)) return(y)
    if (anyDuplicated(sx) || anyDuplicated(sy) || !setequal(sx, sy)) {
        stop("the sample names of `x` and `y` must be unique and identical",
            call. = FALSE)
    }
    y[sx, , drop = FALSE]
}
