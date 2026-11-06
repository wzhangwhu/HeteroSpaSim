# Protocol isolated-v1: stable named streams, independent of traversal order.
# No change to the caller's .Random.seed in isolated mode; legacy evaluates inline.
.with_rng_stream <- function(cfg, key, expr) {
    if (is.null(cfg) || identical(cfg$rng_mode %||% "legacy", "legacy")) return(force(expr))
    had <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
    previous <- if (had) get(".Random.seed", envir = .GlobalEnv) else NULL
    on.exit({
        if (had) assign(".Random.seed", previous, envir = .GlobalEnv)
        else if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    h <- as.double(cfg$seed) %% 2147483646
    for (b in utf8ToInt(paste0("HeteroSpaSim-isolated-v1/", key))) h <- (h * 131 + b) %% 2147483646
    set.seed(as.integer(h + 1))
    force(expr)
}
