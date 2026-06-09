tsrclust <- function(x, k = NULL, alpha = NULL, restr_fact = 20,
                     eps = NULL, min_pts = NULL,
                     macro_method = c("tclust", "tkmeans"),
                     auto_macro = FALSE,
                     k_grid = 1:5,
                     alpha_grid = seq(0, 0.20, length.out = 6),
                     c_grid = c(1, 2, 4, 8, 16, 32, 64, 128),
                     macro_ic = c("MIXMIX", "MIXCLA", "CLACLA"),
                     validated = TRUE,
                     null_reps = 19L,
                     structure_alpha = 0.05,
                     min_density_ratio = 2,
                     novelty_prob = 0.995,
                     eps_quantile = 0.90,
                     standardize = TRUE,
                     ...) {
  call <- match.call()
  x <- .as_numeric_matrix(x)
  macro_method <- match.arg(macro_method)
  macro_ic <- match.arg(macro_ic)
  min_pts <- if (is.null(min_pts)) max(5L, 2L * ncol(x)) else as.integer(min_pts)
  if (min_pts < 2L) {
    stop("min_pts must be at least 2.", call. = FALSE)
  }

  macro <- .fit_macro_stage(
    x = x, k = k, alpha = alpha, restr_fact = restr_fact,
    macro_method = macro_method, auto_macro = auto_macro,
    k_grid = k_grid, alpha_grid = alpha_grid, c_grid = c_grid,
    macro_ic = macro_ic, ...
  )
  macro_labels <- macro$cluster
  main_k <- macro$parameters$k
  final_labels <- as.integer(macro_labels)
  residual_index <- which(final_labels == 0L)
  micro_labels <- integer(length(residual_index))
  micro_info <- list(eps = NA_real_, p_value = NA_real_, score = NA_real_,
                     threshold = NA_real_, n_candidates = 0L)

  if (length(residual_index) >= min_pts) {
    x_distance <- if (standardize) .robust_standardize(x) else x
    residual_x <- x_distance[residual_index, , drop = FALSE]

    if (is.null(eps)) {
      if (validated) {
        macro_reference <- .make_macro_reference(x_distance, macro_labels, main_k)
        novelty_threshold <- stats::qchisq(novelty_prob, df = ncol(x))
        micro_result <- .discover_validated_microclusters(
          residual_x, min_pts = min_pts, null_reps = as.integer(null_reps),
          structure_alpha = structure_alpha,
          min_density_ratio = min_density_ratio,
          macro_reference = macro_reference,
          novelty_threshold = novelty_threshold
        )
      } else {
        selected_eps <- .select_eps_knn(residual_x, min_pts, eps_quantile)
        micro_result <- if (is.na(selected_eps)) {
          list(cluster = integer(nrow(residual_x)), eps = NA_real_)
        } else {
          list(cluster = as.integer(dbscan::dbscan(residual_x, eps = selected_eps, minPts = min_pts)$cluster),
               eps = selected_eps)
        }
      }
    } else {
      micro_result <- list(
        cluster = as.integer(dbscan::dbscan(residual_x, eps = eps, minPts = min_pts)$cluster),
        eps = eps
      )
    }

    micro_labels <- micro_result$cluster
    if (any(micro_labels > 0L)) {
      final_labels[residual_index[micro_labels > 0L]] <- main_k + micro_labels[micro_labels > 0L]
    }
    micro_info <- utils::modifyList(micro_info, micro_result[names(micro_result) != "cluster"])
  }

  out <- list(
    cluster = final_labels,
    raw_cluster = as.integer(macro_labels),
    micro_cluster = micro_labels,
    residual_index = residual_index,
    n_main = main_k,
    n_minority = length(unique(final_labels[final_labels > main_k])),
    n_noise = sum(final_labels == 0L),
    min_pts = min_pts,
    macro_method = macro_method,
    macro_model = macro$model,
    macro_parameters = macro$parameters,
    macro_selection = macro$selection,
    micro_info = micro_info,
    call = call
  )
  class(out) <- "tsrclust"
  out
}

.fit_macro_stage <- function(x, k, alpha, restr_fact, macro_method, auto_macro,
                             k_grid, alpha_grid, c_grid, macro_ic, ...) {
  if (macro_method == "tkmeans" && auto_macro) {
    stop("auto_macro is currently available only for macro_method = 'tclust'.", call. = FALSE)
  }
  if (!auto_macro) {
    if (is.null(k) || is.null(alpha)) {
      stop("k and alpha must be provided when auto_macro = FALSE.", call. = FALSE)
    }
    if (length(k) != 1L || length(alpha) != 1L) {
      stop("k and alpha must be scalar values when auto_macro = FALSE.", call. = FALSE)
    }
    if (macro_method == "tclust") {
      model <- tclust::tclust(x, k = k, alpha = alpha, restr.fact = restr_fact, ...)
      labels <- as.integer(model$cluster)
      parameters <- list(k = as.integer(k), alpha = alpha, restr_fact = restr_fact)
    } else {
      model <- tclust::tkmeans(x, k = k, alpha = alpha, ...)
      labels <- .normalise_tkmeans_labels(model$cluster, as.integer(k))
      parameters <- list(k = as.integer(k), alpha = alpha, restr_fact = NA_real_)
    }
    return(list(model = model, cluster = labels, parameters = parameters, selection = NULL))
  }

  kk <- if (is.null(k)) k_grid else k
  aa <- if (is.null(alpha)) alpha_grid else alpha
  cc <- if (is.null(restr_fact)) c_grid else restr_fact
  candidates <- list()
  for (alpha_i in aa) {
    ic_res <- tclust::tclustIC(x, kk = kk, cc = cc, alpha = alpha_i, whichIC = macro_ic, ...)
    ic_values <- ic_res[[macro_ic]]
    if (is.null(ic_values)) {
      next
    }
    ic_values[!is.finite(ic_values)] <- Inf
    if (all(is.infinite(ic_values))) {
      next
    }
    best_pos <- which(ic_values == min(ic_values), arr.ind = TRUE)[1L, ]
    candidates[[length(candidates) + 1L]] <- list(
      alpha = alpha_i,
      k = ic_res$kk[best_pos[1L]],
      restr_fact = ic_res$cc[best_pos[2L]],
      ic = ic_values[best_pos[1L], best_pos[2L]]
    )
  }
  if (length(candidates) == 0L) {
    stop("No valid TCLUST macro solution was found.", call. = FALSE)
  }
  best <- candidates[[which.min(vapply(candidates, function(z) z$ic, numeric(1L)))]]
  model <- tclust::tclust(x, k = best$k, alpha = best$alpha, restr.fact = best$restr_fact, ...)
  list(
    model = model,
    cluster = as.integer(model$cluster),
    parameters = list(k = as.integer(best$k), alpha = best$alpha, restr_fact = best$restr_fact),
    selection = list(criterion = macro_ic, selected_value = best$ic, candidates = candidates)
  )
}
