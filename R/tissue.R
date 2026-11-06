.simulate_coordinates <- function(n, geometry = "irregular") {
    if (geometry == "irregular") {
        xy <- cbind(stats::runif(n), stats::runif(n))
    } else {
        side <- ceiling(sqrt(n))
        g <- expand.grid(x = seq(0, 1, length.out = side), y = seq(0, 1, length.out = side))
        xy <- as.matrix(g[seq_len(n), c("x", "y")])
        if (geometry == "hex") {
            row_id <- rep(seq_len(side), each = side)[seq_len(n)]
            xy[, 1L] <- xy[, 1L] + (row_id %% 2L) * (0.5 / max(side - 1L, 1L))
            xy[, 1L] <- (xy[, 1L] - min(xy[, 1L])) / max(diff(range(xy[, 1L])), 1e-8)
        }
    }
    colnames(xy) <- c("x", "y")
    xy
}

.simulate_domains <- function(coords, K, pattern, spatial_strength) {
    n <- nrow(coords)
    if (K <= 1L) return(rep.int(1L, n))
    strength <- .clip(spatial_strength, 0, 1)
    noise_sd <- 0.02 + 0.30 * (1 - strength)
    x <- coords[, 1L]
    y <- coords[, 2L]

    if (pattern == "layered") {
        score <- y + stats::rnorm(n, 0, noise_sd)
        br <- stats::quantile(score, probs = seq(0, 1, length.out = K + 1L), names = FALSE)
        br[1L] <- -Inf; br[length(br)] <- Inf
        z <- as.integer(cut(score, breaks = br, include.lowest = TRUE, labels = FALSE))
    } else if (pattern == "nested") {
        r <- sqrt((x - 0.5)^2 + (y - 0.5)^2) + stats::rnorm(n, 0, noise_sd / 2)
        br <- stats::quantile(r, probs = seq(0, 1, length.out = K + 1L), names = FALSE)
        br[1L] <- -Inf; br[length(br)] <- Inf
        z <- as.integer(cut(r, breaks = br, include.lowest = TRUE, labels = FALSE))
    } else if (pattern == "gradient") {
        g <- 0.65 * x + 0.35 * y + stats::rnorm(n, 0, noise_sd)
        br <- stats::quantile(g, probs = seq(0, 1, length.out = K + 1L), names = FALSE)
        br[1L] <- -Inf; br[length(br)] <- Inf
        z <- as.integer(cut(g, breaks = br, include.lowest = TRUE, labels = FALSE))
    } else {
        seed_idx <- sample(seq_len(n), K, replace = FALSE)
        centers <- coords[seed_idx, , drop = FALSE]
        d2 <- sapply(seq_len(K), function(k) rowSums((coords - matrix(centers[k, ], n, 2L, byrow = TRUE))^2))
        dscale <- stats::median(apply(d2, 1L, min)) + 1e-8
        logits <- -d2 / (dscale * (0.35 + 2.5 * (1 - strength)))
        logits <- logits + matrix(stats::rnorm(n * K, 0, 0.25 * (1 - strength)), n, K)
        z <- max.col(logits, ties.method = "random")
    }
    as.integer(z)
}

.simulate_celltypes <- function(domains, K, C, domain_purity = 0.65, composition = NULL) {
    if (!is.null(composition)) {
        comp <- as.matrix(composition)
        comp <- comp / rowSums(comp)
    } else {
        comp <- matrix((1 - domain_purity) / max(C - 1, 1), nrow = K, ncol = C)
        if (C == 1L) comp[,] <- 1
    }
    if (is.null(composition) && C > 1L) {
        for (k in seq_len(K)) {
            pref <- ((k - 1L) %% C) + 1L
            comp[k, pref] <- domain_purity
            comp[k, ] <- 0.85 * comp[k, ] + 0.15 * .rdirichlet1(rep(2, C))
            comp[k, ] <- comp[k, ] / sum(comp[k, ])
        }
    }
    ctype <- integer(length(domains))
    for (k in seq_len(K)) {
        idx <- which(domains == k)
        if (length(idx)) ctype[idx] <- sample.int(C, length(idx), replace = TRUE, prob = comp[k, ])
    }
    list(cell_type = ctype, composition = comp)
}

.assign_batches <- function(coords, B, mode = "random") {
    n <- nrow(coords)
    if (B <= 1L) return(rep.int(1L, n))
    if (mode == "spatial") {
        br <- stats::quantile(coords[, 1L], probs = seq(0, 1, length.out = B + 1L), names = FALSE)
        if (anyDuplicated(br)) return(as.integer(pmin(B, ceiling(rank(coords[, 1L], ties.method = "first") * B / n))))
        br[1L] <- -Inf; br[length(br)] <- Inf
        as.integer(cut(coords[, 1L], breaks = br, include.lowest = TRUE, labels = FALSE))
    } else {
        sample(rep(seq_len(B), length.out = n), n, replace = FALSE)
    }
}

.rbf_basis <- function(coords, n_basis = 25L, length_scale = 0.2) {
    n <- nrow(coords)
    B <- min(as.integer(n_basis), n)
    knot_idx <- unique(round(seq(1, n, length.out = B)))
    if (length(knot_idx) < B) knot_idx <- sample(seq_len(n), B)
    knots <- coords[knot_idx, , drop = FALSE]
    x2 <- rowSums(coords^2)
    k2 <- rowSums(knots^2)
    d2 <- outer(x2, k2, "+") - 2 * coords %*% t(knots)
    d2[d2 < 0] <- 0
    phi <- exp(-d2 / (2 * max(length_scale, 1e-4)^2))
    phi <- sweep(phi, 1L, pmax(rowSums(phi), 1e-10), "/")
    list(phi = phi, knots = knots)
}

.composition_from_labels <- function(domains, celltypes, K = max(domains), C = max(celltypes)) {
    tab <- table(factor(domains, levels = seq_len(K)), factor(celltypes, levels = seq_len(C)))
    comp <- as.matrix(tab) + 1e-8
    comp / rowSums(comp)
}
