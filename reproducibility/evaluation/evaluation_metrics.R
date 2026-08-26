cluster_f1_matrix <- function(truth, predicted, true_labels, pred_labels) {
  if (length(true_labels) == 0L || length(pred_labels) == 0L) {
    return(matrix(0, nrow = length(true_labels), ncol = length(pred_labels),
                  dimnames = list(as.character(true_labels), as.character(pred_labels))))
  }
  out <- matrix(0, nrow = length(true_labels), ncol = length(pred_labels),
                dimnames = list(as.character(true_labels), as.character(pred_labels)))
  for (i in seq_along(true_labels)) {
    ti <- truth == true_labels[i]
    n_t <- sum(ti)
    for (j in seq_along(pred_labels)) {
      pj <- predicted == pred_labels[j]
      overlap <- sum(ti & pj)
      out[i, j] <- safe_ratio(2 * overlap, n_t + sum(pj), 0)
    }
  }
  out
}

max_assignment_dp <- function(score) {
  m <- nrow(score)
  J <- ncol(score)
  if (m == 0L) return(integer(0L))
  if (J == 0L) return(rep(0L, m))
  nstate <- 2^m
  dp <- rep(-Inf, nstate)
  paths <- vector("list", nstate)
  dp[1L] <- 0
  paths[[1L]] <- rep(0L, m)
  for (j in seq_len(J)) {
    new_dp <- dp
    new_paths <- paths
    for (mask in 0:(nstate - 1L)) {
      idx <- mask + 1L
      if (!is.finite(dp[idx])) next
      for (i in seq_len(m)) {
        bit <- bitwShiftL(1L, i - 1L)
        if (bitwAnd(mask, bit) != 0L) next
        newmask <- bitwOr(mask, bit)
        newidx <- newmask + 1L
        candidate <- dp[idx] + score[i, j]
        if (candidate > new_dp[newidx] + 1e-15) {
          new_dp[newidx] <- candidate
          path <- paths[[idx]]
          if (is.null(path)) path <- rep(0L, m)
          path[i] <- j
          new_paths[[newidx]] <- path
        }
      }
    }
    dp <- new_dp
    paths <- new_paths
  }
  best <- which.max(dp)
  path <- paths[[best]]
  if (is.null(path)) rep(0L, m) else path
}

hierarchical_match <- function(truth, predicted, main_k) {
  pred_clusters <- sort(unique(predicted[predicted > 0L]))
  true_major <- seq_len(main_k)
  major_score <- cluster_f1_matrix(truth, predicted, true_major, pred_clusters)
  major_cols <- max_assignment_dp(major_score)
  matched_major_pred <- rep(0L, length(major_cols))
  major_ok <- which(major_cols > 0L)
  if (length(major_ok)) matched_major_pred[major_ok] <- pred_clusters[major_cols[major_ok]]
  used <- matched_major_pred[matched_major_pred > 0L]
  remaining <- setdiff(pred_clusters, used)
  true_minority <- sort(unique(truth[truth > main_k]))
  minority_score <- cluster_f1_matrix(truth, predicted, true_minority, remaining)
  minority_cols <- max_assignment_dp(minority_score)
  matched_minority_pred <- rep(0L, length(minority_cols))
  minority_ok <- which(minority_cols > 0L)
  if (length(minority_ok)) matched_minority_pred[minority_ok] <- remaining[minority_cols[minority_ok]]
  list(
    pred_clusters = pred_clusters,
    true_major = true_major,
    true_minority = true_minority,
    major_score = major_score,
    major_cols = major_cols,
    matched_major_pred = matched_major_pred,
    remaining = remaining,
    minority_score = minority_score,
    minority_cols = minority_cols,
    matched_minority_pred = matched_minority_pred
  )
}

evaluate_partition <- function(truth, predicted, main_k) {
  truth <- as.integer(truth)
  predicted <- as.integer(predicted)
  mt <- hierarchical_match(truth, predicted, main_k)
  major_idx <- truth %in% seq_len(main_k)
  minority_idx <- truth > main_k
  noise_idx <- truth == 0L
  major_f1 <- if (length(mt$true_major) > 0L) {
    mean(vapply(seq_along(mt$true_major), function(i) {
      j <- mt$major_cols[i]
      if (j == 0L) 0 else mt$major_score[i, j]
    }, numeric(1L)))
  } else NA_real_
  minority_cluster_f1 <- if (length(mt$true_minority) > 0L) {
    mean(vapply(seq_along(mt$true_minority), function(i) {
      j <- mt$minority_cols[i]
      if (j == 0L) 0 else mt$minority_score[i, j]
    }, numeric(1L)))
  } else NA_real_
  predicted_positive <- predicted %in% mt$matched_minority_pred[mt$matched_minority_pred > 0L]
  tp <- sum(minority_idx & predicted_positive)
  fp <- sum(!minority_idx & predicted_positive)
  fn <- sum(minority_idx & !predicted_positive)
  tn <- sum(!minority_idx & !predicted_positive)
  precision <- safe_ratio(tp, tp + fp, if (any(minority_idx)) 0 else NA_real_)
  recall <- safe_ratio(tp, tp + fn, if (any(minority_idx)) 0 else NA_real_)
  specificity <- safe_ratio(tn, tn + fp, NA_real_)
  f1 <- if (!any(minority_idx)) NA_real_ else if ((precision + recall) == 0) 0 else 2 * precision * recall / (precision + recall)
  bacc <- if (!any(minority_idx)) NA_real_ else mean(c(recall, specificity), na.rm = TRUE)
  mcc_den <- sqrt(
    as.double(tp + fp) * as.double(tp + fn) *
    as.double(tn + fp) * as.double(tn + fn)
  )
  mcc_num <- as.double(tp) * as.double(tn) -
    as.double(fp) * as.double(fn)
  mcc <- if (!any(minority_idx)) NA_real_ else
    safe_ratio(mcc_num, mcc_den, 0)
  matched_major_union <- mt$matched_major_pred[mt$matched_major_pred > 0L]
  major_retention <- if (any(major_idx)) mean(predicted[major_idx] %in% matched_major_union) else NA_real_
  q_true <- length(mt$true_minority)
  q_hat <- length(mt$remaining)
  matched_minority_nonzero <- mt$matched_minority_pred[mt$matched_minority_pred > 0L]
  spurious <- length(setdiff(mt$remaining, matched_minority_nonzero))
  extra_idx <- predicted %in% mt$remaining
  false_extra_assignment <- if (sum(!minority_idx) > 0L) mean(extra_idx[!minority_idx]) else NA_real_
  minority_ari <- if (q_true >= 2L && sum(minority_idx) >= 2L) adjusted_rand_index(truth[minority_idx], predicted[minority_idx]) else NA_real_
  data.frame(
    GlobalARI = adjusted_rand_index(truth, predicted),
    MajorARI = if (sum(major_idx) >= 2L) adjusted_rand_index(truth[major_idx], predicted[major_idx]) else NA_real_,
    MajorClusterF1 = major_f1,
    MajorRetention = major_retention,
    MinorityPrecision = precision,
    MinorityRecall = recall,
    MinorityDetectionF1 = f1,
    MinorityBACC = bacc,
    MinorityMCC = mcc,
    MinorityClusterF1 = minority_cluster_f1,
    MinorityARI = minority_ari,
    QTrue = q_true,
    QHat = q_hat,
    QError = q_hat - q_true,
    AbsQError = abs(q_hat - q_true),
    SpuriousExtraClusters = spurious,
    AnyMicroDiscovery = as.numeric(q_hat > 0L),
    FalseExtraAssignmentRate = false_extra_assignment,
    NoiseRecall = if (any(noise_idx)) mean(predicted[noise_idx] == 0L) else NA_real_,
    row.names = NULL
  )
}
