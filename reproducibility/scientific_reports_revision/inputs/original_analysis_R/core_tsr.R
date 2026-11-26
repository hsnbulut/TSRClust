select_eps_knn <- function(x, min_pts, eps_quantile = 0.90) {
  if (nrow(x) < min_pts) return(NA_real_)
  k <- max(1L, min(min_pts - 1L, nrow(x) - 1L))
  distances <- dbscan::kNNdist(x, k = k)
  distances <- distances[is.finite(distances) & distances > 0]
  if (length(distances) == 0L) return(NA_real_)
  as.numeric(stats::quantile(distances, probs = eps_quantile, names = FALSE, type = 8))
}

make_eps_grid <- function(x, min_pts, grid_size = 31L) {
  if (nrow(x) < min_pts) return(numeric(0L))
  k <- max(1L, min(min_pts - 1L, nrow(x) - 1L))
  distances <- dbscan::kNNdist(x, k = k)
  distances <- distances[is.finite(distances) & distances > 0]
  if (length(distances) == 0L) return(numeric(0L))
  probs <- seq(0.02, 0.98, length.out = grid_size)
  unique(as.numeric(stats::quantile(distances, probs = probs, names = FALSE, type = 8)))
}

make_macro_reference <- function(x, macro_labels, main_k) {
  lapply(seq_len(main_k), function(cluster_id) {
    cluster_x <- x[macro_labels == cluster_id, , drop = FALSE]
    if (nrow(cluster_x) < 2L) {
      covariance <- diag(ncol(x))
      center <- colMeans(x)
    } else {
      covariance <- stats::cov(cluster_x)
      center <- colMeans(cluster_x)
    }
    if (ncol(x) == 1L) covariance <- matrix(covariance, 1L, 1L)
    ridge <- max(mean(diag(covariance)), 1) * 1e-8
    list(center = center, covariance = covariance + diag(ridge, ncol(x)))
  })
}

candidate_novelty_distance <- function(x, index, macro_reference) {
  if (is.null(macro_reference) || length(macro_reference) == 0L) return(Inf)
  candidate_center <- colMeans(x[index, , drop = FALSE])
  min(vapply(macro_reference, function(reference) {
    stats::mahalanobis(candidate_center, center = reference$center, cov = reference$covariance)
  }, numeric(1L)))
}

cluster_density_scores <- function(x, labels, min_pts, min_density_ratio,
                                   macro_reference, novelty_threshold) {
  cluster_ids <- sort(unique(labels[labels > 0L]))
  if (length(cluster_ids) == 0L) {
    return(data.frame(cluster = integer(0L), size = integer(0L),
                      density_ratio = numeric(0L), novelty_distance = numeric(0L),
                      score = numeric(0L)))
  }
  k <- max(1L, min(min_pts - 1L, nrow(x) - 1L))
  global_distance <- stats::median(dbscan::kNNdist(x, k = k))
  out <- lapply(cluster_ids, function(cluster_id) {
    index <- which(labels == cluster_id)
    if (length(index) < min_pts) {
      return(data.frame(cluster = cluster_id, size = length(index),
                        density_ratio = 0, novelty_distance = 0, score = 0))
    }
    within_k <- max(1L, min(k, length(index) - 1L))
    within_distance <- stats::median(dbscan::kNNdist(x[index, , drop = FALSE], k = within_k))
    density_ratio <- global_distance / max(within_distance, sqrt(.Machine$double.eps))
    novelty_distance <- candidate_novelty_distance(x, index, macro_reference)
    score <- if (density_ratio >= min_density_ratio && novelty_distance >= novelty_threshold) {
      length(index) * log(density_ratio)
    } else 0
    data.frame(cluster = cluster_id, size = length(index),
               density_ratio = density_ratio, novelty_distance = novelty_distance,
               score = score)
  })
  do.call(rbind, out)
}

scan_dbscan_candidates <- function(x, min_pts, eps_grid, min_density_ratio,
                                   macro_reference, novelty_threshold) {
  candidates <- lapply(eps_grid, function(eps) {
    labels <- dbscan::dbscan(x, eps = eps, minPts = min_pts)$cluster
    scores <- cluster_density_scores(x, labels, min_pts, min_density_ratio,
                                     macro_reference, novelty_threshold)
    scores <- scores[scores$score > 0, , drop = FALSE]
    list(eps = eps, labels = labels, scores = scores)
  })
  candidates[vapply(candidates, function(z) nrow(z$scores) > 0L, logical(1L))]
}

max_candidate_score <- function(candidates) {
  if (length(candidates) == 0L) return(0)
  max(vapply(candidates, function(z) max(z$scores$score), numeric(1L)))
}

simulate_diffuse_noise <- function(x) {
  limits <- apply(x, 2L, range)
  out <- vapply(seq_len(ncol(x)), function(j) {
    stats::runif(nrow(x), min = limits[1L, j], max = limits[2L, j])
  }, numeric(nrow(x)))
  matrix(out, nrow = nrow(x), ncol = ncol(x))
}

discover_validated_microclusters <- function(x, min_pts, null_reps,
                                             structure_alpha, min_density_ratio,
                                             macro_reference, novelty_threshold,
                                             eps_grid_size = 31L) {
  empty <- list(cluster = integer(nrow(x)), eps = NA_real_, p_value = 1,
                score = 0, threshold = NA_real_, n_candidates = 0L)
  if (nrow(x) < min_pts) return(empty)
  eps_grid <- make_eps_grid(x, min_pts, eps_grid_size)
  observed <- scan_dbscan_candidates(x, min_pts, eps_grid, min_density_ratio,
                                     macro_reference, novelty_threshold)
  observed_max <- max_candidate_score(observed)
  if (observed_max <= 0 || length(observed) == 0L) return(empty)
  null_max <- numeric(null_reps)
  for (b in seq_len(null_reps)) {
    null_candidates <- scan_dbscan_candidates(simulate_diffuse_noise(x), min_pts, eps_grid,
                                              min_density_ratio, macro_reference,
                                              novelty_threshold)
    null_max[b] <- max_candidate_score(null_candidates)
  }
  p_value <- (1 + sum(null_max >= observed_max)) / (null_reps + 1)
  if (p_value > structure_alpha) {
    empty$p_value <- p_value
    empty$score <- observed_max
    empty$threshold <- max(null_max)
    empty$n_candidates <- sum(vapply(observed, function(z) nrow(z$scores), integer(1L)))
    return(empty)
  }
  threshold <- max(null_max)
  eligible <- lapply(observed, function(z) {
    accepted <- z$scores$score > threshold
    z$accepted_clusters <- z$scores$cluster[accepted]
    z$quality <- sum(z$scores$size[accepted] * z$scores$density_ratio[accepted]^2)
    z$coverage <- sum(z$scores$size[accepted])
    z
  })
  eligible <- eligible[vapply(eligible, function(z) length(z$accepted_clusters) > 0L, logical(1L))]
  if (length(eligible) == 0L) {
    empty$p_value <- p_value
    empty$score <- observed_max
    empty$threshold <- threshold
    return(empty)
  }
  quality <- vapply(eligible, function(z) z$quality, numeric(1L))
  coverage <- vapply(eligible, function(z) z$coverage, numeric(1L))
  best <- eligible[[order(quality, coverage, decreasing = TRUE)[1L]]]
  labels <- best$labels
  labels[!labels %in% best$accepted_clusters] <- 0L
  positive <- sort(unique(labels[labels > 0L]))
  for (j in seq_along(positive)) labels[labels == positive[j]] <- j
  list(cluster = as.integer(labels), eps = best$eps, p_value = p_value,
       score = observed_max, threshold = threshold,
       n_candidates = sum(vapply(observed, function(z) nrow(z$scores), integer(1L))))
}

fit_tsr_from_macro <- function(x, macro_labels, main_k, validated = TRUE,
                               null_reps = 49L, structure_alpha = 0.05,
                               min_density_ratio = 2, novelty_prob = 0.995,
                               eps_quantile = 0.90, eps_grid_size = 31L) {
  final_labels <- as.integer(macro_labels)
  residual_index <- which(final_labels == 0L)
  min_pts <- max(5L, 2L * ncol(x))
  micro <- list(cluster = integer(length(residual_index)), eps = NA_real_,
                p_value = NA_real_, score = NA_real_, threshold = NA_real_,
                n_candidates = 0L)
  if (length(residual_index) >= min_pts) {
    residual_x <- x[residual_index, , drop = FALSE]
    if (validated) {
      macro_reference <- make_macro_reference(x, macro_labels, main_k)
      novelty_threshold <- stats::qchisq(novelty_prob, df = ncol(x))
      micro <- discover_validated_microclusters(
        residual_x, min_pts, as.integer(null_reps), structure_alpha,
        min_density_ratio, macro_reference, novelty_threshold, eps_grid_size
      )
    } else {
      eps <- select_eps_knn(residual_x, min_pts, eps_quantile)
      if (is.finite(eps)) {
        micro <- list(cluster = as.integer(dbscan::dbscan(residual_x, eps = eps,
                                                         minPts = min_pts)$cluster),
                      eps = eps, p_value = NA_real_, score = NA_real_,
                      threshold = NA_real_, n_candidates = NA_integer_)
      }
    }
    positive <- micro$cluster > 0L
    if (any(positive)) final_labels[residual_index[positive]] <- main_k + micro$cluster[positive]
  }
  list(cluster = final_labels, raw_cluster = as.integer(macro_labels),
       residual_index = residual_index, micro = micro, min_pts = min_pts,
       n_minority = length(unique(final_labels[final_labels > main_k])))
}
