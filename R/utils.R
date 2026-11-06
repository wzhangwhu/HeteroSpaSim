`%||%` <- function(x, y) {
    if (is.null(x)) y else x
}

.clip <- function(x, lower, upper) {
    pmin(pmax(x, lower), upper)
}

.recycle_named <- function(x, names_target, default = NULL) {
    if (is.null(x)) x <- default
    if (length(x) == 1L) {
        out <- rep(x, length(names_target))
        names(out) <- names_target
        return(out)
    }
    if (is.null(names(x))) {
        if (length(x) != length(names_target)) {
            stop("Unnamed parameter length must match the number of modalities.", call. = FALSE)
        }
        names(x) <- names_target
        return(x)
    }
    missing_names <- setdiff(names_target, names(x))
    if (length(missing_names)) {
        stop("Missing values for modalities: ", paste(missing_names, collapse = ", "), call. = FALSE)
    }
    x[names_target]
}

.standardize_cols <- function(x, eps = 1e-8) {
    x <- as.matrix(x)
    mu <- colMeans(x)
    x <- sweep(x, 2L, mu, "-")
    ss <- sqrt(colMeans(x^2))
    ss[!is.finite(ss) | ss < eps] <- 1
    sweep(x, 2L, ss, "/")
}

.rdirichlet1 <- function(alpha) {
    x <- stats::rgamma(length(alpha), shape = alpha, rate = 1)
    if (sum(x) <= 0) rep(1 / length(alpha), length(alpha)) else x / sum(x)
}

.row_mode_from_operator <- function(A, labels) {
    if (!inherits(A, "Matrix")) A <- Matrix::Matrix(A, sparse = TRUE)
    lev <- sort(unique(labels))
    out <- integer(nrow(A))
    for (j in seq_len(nrow(A))) {
        idx <- which(A[j, ] != 0)
        if (!length(idx)) {
            out[j] <- NA_integer_
        } else {
            labs <- labels[idx]
            ww <- as.numeric(A[j, idx])
            tab <- vapply(lev, function(v) sum(ww[labs == v]), numeric(1))
            out[j] <- lev[which.max(tab)]
        }
    }
    out
}

.sparse_cbind <- function(blocks) {
    if (!length(blocks)) return(Matrix::Matrix(nrow = 0, ncol = 0, sparse = TRUE))
    if (length(blocks) == 1L) return(blocks[[1L]])
    Reduce(function(a, b) cbind(a, b), blocks)
}

.sample_without <- function(pool, n, exclude = integer()) {
    pool <- setdiff(pool, exclude)
    if (n <= 0L) return(integer())
    if (length(pool) < n) n <- length(pool)
    if (n == 0L) return(integer())
    pool[sample.int(length(pool), n, replace = FALSE)]
}

.safe_dir <- function(path) {
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    normalizePath(path, winslash = "/", mustWork = TRUE)
}

.normalize_coordinates <- function(xy) {
    xy <- as.matrix(xy)
    if (ncol(xy) < 2L) stop("coordinates must have at least two columns.", call. = FALSE)
    xy <- xy[, 1:2, drop = FALSE]
    storage.mode(xy) <- "double"
    if (any(!is.finite(xy))) stop("coordinates must be finite.", call. = FALSE)
    for (j in 1:2) {
        r <- range(xy[, j], finite = TRUE)
        if (!all(is.finite(r))) stop("coordinates contain non-finite values.", call. = FALSE)
        d <- diff(r)
        xy[, j] <- if (d < 1e-12) 0.5 else (xy[, j] - r[1L]) / d
    }
    colnames(xy) <- c("x", "y")
    xy
}
