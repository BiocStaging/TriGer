MATLAB_GREY <- grDevices::rgb(0.15, 0.15, 0.15)

#' Jet colour map
#'
#' The jet colour map used by the plotting functions.
#'
#' @param m Number of colours.
#' @return A character vector of `m` colours.
#' @examples
#' matlab_jet(8)
#' @export
matlab_jet <- function(m = 256) {
  n <- ceiling(m / 4)
  u <- c((1:n) / n, rep(1, n - 1), (n:1) / n)
  g <- ceiling(n / 2) - (m %% 4 == 1) + seq_along(u)
  r <- g + n; b <- g - n
  r <- r[r <= m]; g <- g[g <= m]; b <- b[b >= 1]
  J <- matrix(0, m, 3)
  J[r, 1] <- u[seq_along(r)]
  J[g, 2] <- u[seq_along(g)]
  J[b, 3] <- u[(length(u) - length(b) + 1):length(u)]
  rgb(J[, 1], J[, 2], J[, 3])
}

# MATLAB CDataMapping 'scaled': index = fix((c - cmin)/(cmax - cmin) * m) + 1, clamped
map_col <- function(v, clim, cmap) {
  m <- length(cmap)
  idx <- floor((v - clim[1]) / (clim[2] - clim[1]) * m) + 1
  cmap[pmin(pmax(idx, 1), m)]
}

unique_keep_last <- function(v) rev(rev(v)[!duplicated(rev(v))])

# Draw many 4-vertex polygons in one call (NA-separated), border = fill to hide seams
draw_quads <- function(px, py, cols) {       # px, py: 4 x k matrices
  k <- ncol(px)
  polygon(as.vector(rbind(px, NA))[-(5 * k)], as.vector(rbind(py, NA))[-(5 * k)],
          col = cols, border = cols, lwd = 0.25)
}

#' Plot detected modules with their within-omics correlations
#'
#' Draws the cross-omics matrix restricted to the selected features as a
#' heatmap, with the X-X and Y-Y matrices of the same features as rotated
#' triangles on the top and right. Draws into the current graphics device;
#' a 5.6 x 4.2 inch device matches the default layout.
#'
#' @param X,Y Within-omics matrices for the X and Y features.
#' @param XY Cross-omics matrix (X features in rows).
#' @param result A list of blocks, each with X indices first and Y indices
#'   second, as returned by [greedy_3step()] or [greedy_peeling_XY_all()].
#' @param line Draw a box around each block.
#' @param step `c(x, y)` tick spacing; defaults to `c(100, 10)`.
#' @param bar Add a colour bar for the middle panel.
#' @param figsize Device size in inches, used to keep the panels square.
#' @param cbar_pos Colour bar position `c(left, bottom, width, height)` in
#'   figure units, or `"right"` to place it beside the plot.
#' @param clim_tri,clim_mid Colour limits for the triangles and the middle
#'   panel. By default `c(-1, 1)` for the triangles, and `c(0, 10)` for the
#'   middle panel if its maximum exceeds 1, otherwise `c(-1, 1)`.
#' @param box_lw Line width of the block boxes, in points.
#' @return Invisibly, a list with the plotted X and Y indices.
#' @examples
#' set.seed(1)
#' N <- 150; m <- 80; n <- 60
#' f <- rnorm(N)
#' X <- matrix(rnorm(N * m), N); Y <- matrix(rnorm(N * n), N)
#' X[, 1:25] <- X[, 1:25] + 1.5 * f; Y[, 1:20] <- Y[, 1:20] + 1.5 * f
#' R <- cor(cbind(X, Y)); P <- abs(R); diag(P) <- 0
#' modules <- greedy_3step(P, R, m, n, 1.1, 1.8, 1.8)
#' plot3in1_v1(R[1:m, 1:m], R[m + 1:n, m + 1:n], R[1:m, m + 1:n], modules,
#'             line = TRUE, step = c(10, 5))
#' @export
plot3in1_v1 <- function(X, Y, XY, result, line = FALSE, step = NULL, bar = FALSE,
                        figsize = c(5.6, 4.2), cbar_pos = c(0.15, 0.13, 0.013, 0.45),
                        clim_tri = NULL, clim_mid = NULL, box_lw = 2) {
  result <- lapply(result, function(r) list(as.integer(unlist(r[[1]])), as.integer(unlist(r[[2]]))))
  select_X <- unique_keep_last(unlist(lapply(result, `[[`, 1)))
  select_Y <- unique_keep_last(unlist(lapply(result, `[[`, 2)))

  X <- as.matrix(X); Y <- as.matrix(Y); XY <- as.matrix(XY)
  C <- Y[select_Y, select_Y, drop = FALSE]
  B <- X[select_X, select_X, drop = FALSE]
  A <- XY[select_X, select_Y, drop = FALSE]
  m <- nrow(A); n <- ncol(A)
  ratio <- m / n / 1.2
  cmap <- matlab_jet(256)
  clim1 <- if (is.null(clim_tri)) c(-1, 1) else clim_tri
  clim2 <- if (!is.null(clim_mid)) clim_mid else if (max(A, na.rm = TRUE) > 1) c(0, 10) else c(-1, 1)

  # ---- axis equal tight: shrink the plot box inside the MATLAB axes Position ----
  xlim <- c(0, max(2 * m / ratio, n + 2 * m / ratio))
  ylim <- c(0, 2 * n + m / ratio)
  pos <- c(0.13, 0.11, 0.775, 0.815)
  box_w <- pos[3] * figsize[1]; box_h <- pos[4] * figsize[2]
  s <- min(box_w / diff(xlim), box_h / diff(ylim))         # inches per data unit
  w <- s * diff(xlim) / figsize[1]; h <- s * diff(ylim) / figsize[2]
  cx <- pos[1] + pos[3] / 2; cy <- pos[2] + pos[4] / 2
  op <- par(plt = c(cx - w / 2, cx + w / 2, cy - h / 2, cy + h / 2), xaxs = "i", yaxs = "i",
            ps = 10, family = "sans", col.axis = MATLAB_GREY, fg = MATLAB_GREY)
  on.exit(par(op), add = TRUE)
  plot.new(); plot.window(xlim, ylim)

  # ---- Clean tick labels (inward ticks, no axis line, MATLAB TickLength 0.01) ----
  # drawn BEFORE the patches: MATLAB Axes Layer = 'bottom' puts tick marks under the data
  xs <- if (is.null(step)) 100 else step[1]; ys <- if (is.null(step)) 10 else step[2]
  xTickIdx <- seq(xs, floor(m / xs) * xs, by = xs)
  yTickIdx <- rev(seq(ys, floor(n / ys) * ys, by = ys))
  tick_in <- 0.01 * max(w * figsize[1], h * figsize[2])     # inches
  tcl <- tick_in / par("csi")
  axis(1, at = (xTickIdx - 1) * 2 / ratio, labels = xTickIdx, lwd = 0, lwd.ticks = 0.5 / 0.75,
       tcl = tcl, mgp = c(0, -0.11, 0), col.ticks = MATLAB_GREY)
  axis(2, at = (n - yTickIdx) * 2 + 1, labels = yTickIdx, lwd = 0, lwd.ticks = 0.5 / 0.75,
       tcl = tcl, mgp = c(3, 0.23, 0), las = 1, col.ticks = MATLAB_GREY)

  # ---- Right: diamonds for Y ----
  ij <- do.call(rbind, lapply(1:n, function(i) cbind(i, i:n)))
  x0 <- ij[, 1] - 1 + m * 2 / ratio; y0 <- 2 * n - ij[, 2] * 2 + ij[, 1]
  vx <- c(0, 1, 0, -1); vy <- c(1, 0, -1, 0)
  draw_quads(outer(vx, x0, `+`), outer(vy, y0, `+`),
             map_col(C[cbind(ij[, 2], ij[, 2] - ij[, 1] + 1)], clim1, cmap))

  # ---- Top: diamonds for X ----
  ij <- do.call(rbind, lapply(1:m, function(i) cbind(i, i:m)))
  x0 <- (ij[, 2] * 2 - ij[, 1]) / ratio; y0 <- (ij[, 1] - 1) / ratio + n * 2
  draw_quads(outer(vx / ratio, x0, `+`), outer(vy / ratio, y0, `+`),
             map_col(B[cbind(ij[, 2], ij[, 2] - ij[, 1] + 1)], clim1, cmap))

  # ---- Middle square (drawn last = on top, as ax2 in MATLAB) ----
  ij <- expand.grid(j = 1:n, i = 1:m)
  x0 <- (ij$i - 1) * 2 / ratio; y0 <- (ij$j - 1) * 2
  cols <- map_col(A[cbind(ij$i, n + 1 - ij$j)], clim2, cmap)
  rect(x0, y0, x0 + 2 / ratio, y0 + 2, col = cols, border = cols, lwd = 0.25)

  # ---- Lines separating blocks ----
  if (line) {
    L <- t(sapply(result, function(r) c(length(r[[1]]), length(r[[2]]))))
    L[, 1] <- L[, 1] * 2 / ratio; L[, 2] <- L[, 2] * 2
    x <- 0; y <- sum(L[, 2]); i <- 1
    while (i <= nrow(L) && L[i, 1] != 0) {
      rect(x, y - L[i, 2], x + L[i, 1], y, border = "red", lwd = box_lw / 0.75)
      y <- y - L[i, 2]; x <- x + L[i, 1]; i <- i + 1
    }
  }

  # ---- Colorbar (middle panel), fixed figure position like MATLAB cb.Position ----
  if (bar) {
    if (is.character(cbar_pos) && cbar_pos == "right") {   # right of the plot, aligned with square
      sq_h <- h * (2 * n) / (2 * n + m / ratio)
      cbar_pos <- c(cx + w / 2 + 0.02, cy - h / 2, 0.013, sq_h)
    }
    par(plt = c(cbar_pos[1], cbar_pos[1] + cbar_pos[3], cbar_pos[2], cbar_pos[2] + cbar_pos[4]),
        new = TRUE)
    plot.new(); plot.window(c(0, 1), clim2)
    yb <- seq(clim2[1], clim2[2], length.out = length(cmap) + 1)
    rect(0, yb[-length(yb)], 1, yb[-1], col = cmap, border = NA)
    axis(4, at = pretty(clim2), labels = format(pretty(clim2), drop0trailing = TRUE, trim = TRUE),
         lwd = 0, lwd.ticks = 0.5 / 0.75, tcl = tcl, las = 1, mgp = c(3, 0.28, 0),
         col.ticks = MATLAB_GREY)
    box(lwd = 0.5 / 0.75, col = MATLAB_GREY)
  }
  invisible(list(select_X = select_X, select_Y = select_Y))
}


#' Heatmap of detected blocks
#'
#' Shows `Wp` reordered so that the features of the detected blocks come
#' first. With `full = TRUE` the whole matrix is shown with one box around the
#' selected region; otherwise only the selected region is shown, with one box
#' per block.
#'
#' @param Wp The matrix the blocks were detected in.
#' @param result A list of blocks, each with X indices first and Y indices
#'   second; for a symmetric `Wp`, only the first element is used.
#' @param full Show the full matrix.
#' @param line Draw boxes around the blocks.
#' @param crange Colour limits; defaults to the range of the plotted values.
#' @return Invisibly, the plotted matrix.
#' @examples
#' set.seed(1)
#' W <- matrix(runif(60 * 80, 0, 0.5), 60, 80)
#' W[6:20, 11:30] <- W[6:20, 11:30] + 4
#' plot_all_result(W, greedy_peeling_XY_all(W, 1.5), full = TRUE)
#' @export
plot_all_result <- function(Wp, result, full = FALSE, line = TRUE, crange = NULL) {
  Wp <- as.matrix(Wp); m <- nrow(Wp); n <- ncol(Wp)
  rows <- lapply(result, function(r) as.integer(r[[1]]))
  cols <- if (m == n) rows else lapply(result, function(r) as.integer(r[[2]]))
  keep <- lengths(rows) > 0; rows <- rows[keep]; cols <- cols[keep]
  sX <- unlist(rows); sY <- unlist(cols)
  M <- if (full) Wp[c(sX, setdiff(seq_len(m), sX)), c(sY, setdiff(seq_len(n), sY)), drop = FALSE]
       else Wp[sX, sY, drop = FALSE]
  if (is.null(crange)) crange <- range(M, na.rm = TRUE)
  M <- pmin(pmax(M, crange[1]), crange[2])
  H <- nrow(M); cmap <- matlab_jet(256)

  op <- par(no.readonly = TRUE); on.exit(par(op))
  layout(matrix(1:2, 1), widths = c(1, graphics::lcm(2.2)))
  par(mar = c(3, 3.5, 1, 0.5))
  image(seq_len(ncol(M)), seq_len(H), t(M[H:1, , drop = FALSE]), col = cmap, zlim = crange,
        axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
  axis(1)
  yt <- pretty(c(1, H)); yt <- yt[yt >= 1 & yt <= H]
  axis(2, at = H + 1 - yt, labels = yt, las = 1)
  if (line) {
    if (full) {
      rect(0.5, H - length(sX) + 0.5, length(sY) + 0.5, H + 0.5, border = "red", lwd = 2)
    } else {
      r0 <- 0; c0 <- 0
      for (k in seq_along(rows)) {
        a <- length(rows[[k]]); b <- length(cols[[k]])
        rect(c0 + 0.5, H - r0 - a + 0.5, c0 + b + 0.5, H - r0 + 0.5, border = "red", lwd = 2)
        r0 <- r0 + a; c0 <- c0 + b
      }
    }
  }
  par(mar = c(3, 0.3, 1, 2.5))
  z <- seq(crange[1], crange[2], length.out = 256)
  image(1, z, matrix(z, 1), col = cmap, axes = FALSE, xlab = "", ylab = "")
  axis(4, las = 1); box()
  invisible(M)
}
