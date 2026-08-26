fit_tclust_fixed <- function(x, k, alpha, restr_fact, nstart, cfg) {
  tclust::tclust(
    x, k = as.integer(k), alpha = alpha, restr.fact = restr_fact,
    nstart = as.integer(nstart), niter1 = cfg$tclust_niter1,
    niter2 = cfg$tclust_niter2, nkeep = cfg$tclust_nkeep,
    opt = "HARD", equal.weights = FALSE, parallel = FALSE,
    store_x = FALSE, trace = 0
  )
}

fit_tkmeans_fixed <- function(x, k, alpha, nstart, cfg) {
  tclust::tkmeans(
    x, k = as.integer(k), alpha = alpha,
    nstart = as.integer(nstart), niter1 = cfg$tclust_niter1,
    niter2 = cfg$tclust_niter2, nkeep = cfg$tclust_nkeep,
    parallel = FALSE, store_x = FALSE, trace = 0
  )
}

fit_auto_tclust <- function(x, cfg, nstart) {
  candidates <- list()
  for (alpha_i in cfg$auto_alpha_grid) {
    ic <- tryCatch(
      tclust::tclustIC(
        x, kk = cfg$auto_k_grid, cc = cfg$auto_c_grid,
        alpha = alpha_i, whichIC = cfg$auto_ic,
        nstart = as.integer(nstart), niter1 = cfg$tclust_niter1,
        niter2 = cfg$tclust_niter2, nkeep = cfg$tclust_nkeep,
        parallel = FALSE, trace = FALSE, store_x = FALSE
      ),
      error = function(e) NULL
    )
    if (is.null(ic) || is.null(ic[[cfg$auto_ic]])) next
    values <- ic[[cfg$auto_ic]]
    values[!is.finite(values)] <- Inf
    if (all(is.infinite(values))) next
    pos <- which(values == min(values), arr.ind = TRUE)[1L, ]
    candidates[[length(candidates) + 1L]] <- list(
      alpha = alpha_i,
      k = ic$kk[pos[1L]],
      restr_fact = ic$cc[pos[2L]],
      ic = values[pos[1L], pos[2L]]
    )
  }
  if (length(candidates) == 0L) stop("No valid automatic TCLUST solution.")
  best <- candidates[[which.min(vapply(candidates, function(z) z$ic, numeric(1L)))]]
  model <- fit_tclust_fixed(x, best$k, best$alpha, best$restr_fact, nstart, cfg)
  list(model = model, selected = best, candidates = candidates)
}

fit_full_dbscan <- function(x, eps_quantile = 0.90) {
  min_pts <- max(5L, 2L * ncol(x))
  eps <- select_eps_knn(x, min_pts, eps_quantile)
  if (!is.finite(eps)) return(list(cluster = integer(nrow(x)), eps = NA_real_, min_pts = min_pts))
  fit <- dbscan::dbscan(x, eps = eps, minPts = min_pts)
  list(cluster = as.integer(fit$cluster), eps = eps, min_pts = min_pts)
}

fit_full_hdbscan <- function(x) {
  min_pts <- max(5L, 2L * ncol(x))
  fit <- dbscan::hdbscan(x, minPts = min_pts)
  list(cluster = as.integer(fit$cluster), min_pts = min_pts)
}
