# Internal input checks shared by the exported functions.

.as_numeric_matrix <- function(x, arg) {
    x <- as.matrix(x)
    if (!is.numeric(x) && !is.logical(x)) {
        stop("`", arg, "` must be a numeric matrix", call. = FALSE)
    }
    storage.mode(x) <- "double"
    if (nrow(x) == 0L || ncol(x) == 0L) {
        stop("`", arg, "` must have at least one row and one column",
            call. = FALSE)
    }
    if (any(!is.finite(x))) {
        stop("`", arg, "` must not contain missing or infinite values",
            call. = FALSE)
    }
    x
}

.check_symmetric <- function(x, arg) {
    if (nrow(x) != ncol(x) ||
        !isTRUE(all.equal(x, t(x), tolerance = 1e-8,
            check.attributes = FALSE))) {
        stop("`", arg, "` must be symmetric", call. = FALSE)
    }
    invisible(x)
}

.check_number <- function(x, arg) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x)) {
        stop("`", arg, "` must be a finite number", call. = FALSE)
    }
    invisible(x)
}

.check_count <- function(x, arg, min = 1) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        x != round(x) || x < min) {
        stop("`", arg, "` must be a whole number of at least ", min,
            call. = FALSE)
    }
    invisible(x)
}

.check_size <- function(x, arg) {
    if (!is.numeric(x) || length(x) != 2L || any(!is.finite(x)) ||
        x[1] > x[2]) {
        stop("`", arg, "` must be c(min, max) with min <= max",
            call. = FALSE)
    }
    invisible(x)
}

.feature_names <- function(x, prefix) {
    nm <- colnames(x)
    if (is.null(nm)) nm <- paste0(prefix, seq_len(ncol(x)))
    nm
}
