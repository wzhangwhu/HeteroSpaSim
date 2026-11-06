.choose2 <- function(x) x * (x - 1) / 2

#' Score spatial-domain recovery
#'
#' @param truth Ground-truth labels.
#' @param predicted Predicted labels.
#' @return Named numeric vector containing ARI, NMI and Purity.
#' @export
score_domains <- function(truth, predicted) {
    if (length(truth) != length(predicted)) stop("truth and predicted must have the same length.", call. = FALSE)
    ok <- !is.na(truth) & !is.na(predicted)
    truth <- as.factor(truth[ok]); predicted <- as.factor(predicted[ok])
    if (length(truth) != length(predicted) || length(truth) < 2L) stop("Need at least two paired labels.", call. = FALSE)
    tab <- table(truth, predicted)
    n <- sum(tab)
    sum_nij <- sum(.choose2(tab))
    sum_ai <- sum(.choose2(rowSums(tab)))
    sum_bj <- sum(.choose2(colSums(tab)))
    total <- .choose2(n)
    expected <- if (total > 0) sum_ai * sum_bj / total else 0
    max_index <- 0.5 * (sum_ai + sum_bj)
    ari <- if (abs(max_index - expected) < 1e-12) 1 else (sum_nij - expected) / (max_index - expected)

    p <- tab / n
    pi <- rowSums(p); pj <- colSums(p)
    nz <- which(p > 0, arr.ind = TRUE)
    mi <- sum(vapply(seq_len(nrow(nz)), function(k) {
        i <- nz[k, 1L]; j <- nz[k, 2L]
        p[i, j] * log(p[i, j] / (pi[i] * pj[j]))
    }, numeric(1)))
    hi <- -sum(pi[pi > 0] * log(pi[pi > 0]))
    hj <- -sum(pj[pj > 0] * log(pj[pj > 0]))
    nmi <- if (hi <= 1e-12 && hj <= 1e-12) {
        1
    } else if (hi <= 1e-12 || hj <= 1e-12) {
        0
    } else {
        mi / sqrt(hi * hj)
    }
    purity <- sum(apply(tab, 2L, max)) / n
    c(ARI = unname(ari), NMI = unname(nmi), Purity = unname(purity))
}

#' Score spatial alignment
#'
#' @param truth True n-by-2 coordinates.
#' @param estimated Estimated n-by-2 coordinates in the same row order.
#' @return Named numeric vector of RMSE, median Euclidean error and 95th percentile error.
#' @export
score_alignment <- function(truth, estimated) {
    truth <- as.matrix(truth); estimated <- as.matrix(estimated)
    if (ncol(truth) < 2L || ncol(estimated) < 2L || !identical(dim(truth), dim(estimated))) stop("truth and estimated coordinates must have matching dimensions with at least two columns.", call. = FALSE)
    if (any(!is.finite(truth[, 1:2, drop = FALSE])) || any(!is.finite(estimated[, 1:2, drop = FALSE]))) stop("Coordinates must be finite.", call. = FALSE)
    d <- sqrt(rowSums((truth[, 1:2, drop = FALSE] - estimated[, 1:2, drop = FALSE])^2))
    c(RMSE = sqrt(mean(d^2)), MedianError = stats::median(d), Q95Error = unname(stats::quantile(d, 0.95)))
}

#' Score continuous-valued imputation
#'
#' @param truth Ground-truth numeric matrix or vector.
#' @param predicted Predicted values with matching dimensions.
#' @return RMSE, Pearson correlation and Spearman correlation.
#' @export
score_imputation <- function(truth, predicted) {
    if (!all(dim(as.matrix(truth)) == dim(as.matrix(predicted)))) stop("truth and predicted must have matching dimensions.", call. = FALSE)
    a <- as.numeric(truth); b <- as.numeric(predicted)
    ok <- is.finite(a) & is.finite(b)
    a <- a[ok]; b <- b[ok]
    if (length(a) < 2L) stop("Not enough finite paired values.", call. = FALSE)
    pear <- suppressWarnings(stats::cor(a, b, method = "pearson"))
    spear <- suppressWarnings(stats::cor(a, b, method = "spearman"))
    c(RMSE = sqrt(mean((a - b)^2)), Pearson = pear, Spearman = spear)
}

#' Score binary recovery such as ATAC accessibility
#'
#' @param truth Binary 0/1 truth.
#' @param score Predicted probabilities or continuous scores.
#' @return AUROC and average-precision AUPRC.
#' @export
score_binary_imputation <- function(truth, score) {
    if (length(truth) != length(score)) stop("truth and score must have the same number of elements.", call. = FALSE)
    if (any(!is.na(truth) & (!is.finite(truth) | !truth %in% c(0, 1)))) stop("truth must be binary 0/1.", call. = FALSE)
    y <- as.integer(truth)
    s <- as.numeric(score)
    ok <- is.finite(y) & is.finite(s)
    y <- y[ok]; s <- s[ok]
    npos <- sum(y == 1L); nneg <- sum(y == 0L)
    if (npos == 0L || nneg == 0L) stop("Both positive and negative classes are required.", call. = FALSE)
    r <- rank(s, ties.method = "average")
    auc <- (sum(r[y == 1L]) - npos * (npos + 1) / 2) / (npos * nneg)
    ord <- order(s, decreasing = TRUE)
    yy <- y[ord]
    tp <- cumsum(yy == 1L)
    fp <- cumsum(yy == 0L)
    precision <- tp / pmax(tp + fp, 1)
    # Evaluate once per distinct score threshold; ties are a single prediction set.
    ends <- c(which(diff(s[ord]) != 0), length(s))
    recall <- tp[ends] / npos
    auprc <- sum(diff(c(0, recall)) * precision[ends])
    c(AUROC = auc, AUPRC = auprc)
}

#' Score cross-modal correspondence recovery
#'
#' @param truth A sparse or dense source-by-target correspondence matrix, where larger values indicate greater true overlap.
#' @param predicted Integer vector giving one predicted target index for each source row.
#' @return Exact top-1 accuracy and mean true overlap at the predicted match.
#' @export
score_correspondence <- function(truth, predicted) {
    C <- truth
    nr <- nrow(C); nc <- ncol(C)
    if (length(predicted) != nr) stop("predicted must contain one target index per source row.", call. = FALSE)
    predicted <- as.integer(predicted)
    if (any(predicted < 1L | predicted > nc, na.rm = TRUE)) stop("predicted contains target indices outside the truth matrix.", call. = FALSE)
    max_truth <- numeric(nr); ov <- numeric(nr)
    for (i in seq_len(nr)) {
        rowi <- as.numeric(C[i, ])
        max_truth[i] <- if (length(rowi)) max(rowi) else 0
        if (!is.na(predicted[i])) ov[i] <- rowi[predicted[i]] else ov[i] <- NA_real_
    }
    ok <- !is.na(predicted) & max_truth > 0
    if (!any(ok)) return(c(Top1Accuracy = NA_real_, MeanTrueOverlap = NA_real_))
    # Any target attaining the row-wise maximum is a valid top-1 match. This
    # avoids penalizing ties under coarse/cross-resolution measurement designs.
    correct <- ov[ok] >= max_truth[ok] - 1e-12
    c(Top1Accuracy = mean(correct), MeanTrueOverlap = mean(ov[ok], na.rm = TRUE))
}

#' Score latent-factor subspace recovery
#'
#' @param truth Ground-truth n-by-q latent representation.
#' @param estimated Estimated n-by-r representation for the same units.
#' @return Mean canonical subspace correlation and normalized projection distance.
#' @export
score_latent_recovery <- function(truth, estimated) {
    X <- as.matrix(truth)
    Y <- as.matrix(estimated)
    if (nrow(X) != nrow(Y)) stop("truth and estimated must have the same number of rows.", call. = FALSE)
    if (ncol(X) == 0L || ncol(Y) == 0L || any(!is.finite(X)) || any(!is.finite(Y))) {
        stop("Latent matrices must be non-empty and finite.", call. = FALSE)
    }

    # Subspace recovery is invariant to translation, so center columns first.
    X <- sweep(X, 2L, colMeans(X), "-")
    Y <- sweep(Y, 2L, colMeans(Y), "-")
    if (qr(X)$rank == 0L || qr(Y)$rank == 0L || !any(abs(X) > 1e-12) || !any(abs(Y) > 1e-12)) stop("Latent matrices must have non-zero rank after centering.", call. = FALSE)

    # Exact self-comparison is a fundamental identity and should not acquire
    # an artificial ~1e-8 distance from QR/SVD round-off.
    if (identical(dim(X), dim(Y)) &&
        isTRUE(all.equal(X, Y, tolerance = 100 * .Machine$double.eps,
                         check.attributes = FALSE))) {
        return(c(MeanCanonicalCorrelation = 1,
                 NormalizedProjectionDistance = 0))
    }

    orth_basis <- function(Z) {
        q <- qr(Z, tol = sqrt(.Machine$double.eps), LAPACK = FALSE)
        if (q$rank <= 0L) return(matrix(numeric(), nrow(Z), 0L))
        qr.Q(q, complete = FALSE)[, seq_len(q$rank), drop = FALSE]
    }

    Qx <- orth_basis(X)
    Qy <- orth_basis(Y)
    if (ncol(Qx) == 0L || ncol(Qy) == 0L) {
        stop("Latent matrices must have non-zero rank after centering.", call. = FALSE)
    }

    sv <- svd(crossprod(Qx, Qy), nu = 0, nv = 0)$d
    sv <- pmin(pmax(sv, 0), 1)

    # Principal correlations that differ from 0 or 1 only at floating-point
    # precision should be treated as exact boundary values. This stabilizes
    # identical and rotation-equivalent subspaces across BLAS/LAPACK builds.
    tol_sv <- max(1e-12, 100 * .Machine$double.eps *
                      max(ncol(Qx), ncol(Qy), nrow(Qx), 1L))
    sv[abs(sv - 1) <= tol_sv] <- 1
    sv[abs(sv) <= tol_sv] <- 0

    rx <- ncol(Qx)
    ry <- ncol(Qy)
    proj_sq <- rx + ry - 2 * sum(sv^2)
    # Small negative values can arise from rounding after summation.
    proj_sq <- max(proj_sq, 0)
    norm_dist <- sqrt(proj_sq / max(rx + ry, 1L))
    if (norm_dist <= sqrt(.Machine$double.eps)) norm_dist <- 0

    c(MeanCanonicalCorrelation = mean(sv),
      NormalizedProjectionDistance = norm_dist)
}
