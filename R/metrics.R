evaluate_tsrclust <- function(truth, predicted, main_k = 2,
                              predicted_positive = predicted > main_k) {
  truth <- as.integer(truth)
  predicted <- as.integer(predicted)
  if (length(truth) != length(predicted)) {
    stop("truth and predicted must have the same length.", call. = FALSE)
  }
  main_index <- truth %in% seq_len(main_k)
  minority_index <- truth > main_k
  noise_index <- truth == 0L
  actual_positive <- minority_index

  tp <- sum(actual_positive & predicted_positive)
  fp <- sum(!actual_positive & predicted_positive)
  fn <- sum(actual_positive & !predicted_positive)
  tn <- sum(!actual_positive & !predicted_positive)
  precision <- .safe_ratio(tp, tp + fp, undefined = if (sum(actual_positive) > 0L) 0 else NA_real_)
  recall <- .safe_ratio(tp, tp + fn)
  specificity <- .safe_ratio(tn, tn + fp)
  f1 <- if (is.na(recall)) {
    NA_real_
  } else if (is.na(precision) || (precision + recall) == 0) {
    0
  } else {
    2 * precision * recall / (precision + recall)
  }

  data.frame(
    GlobalARI = .adjusted_rand_index(truth, predicted),
    MainARI = if (any(main_index)) .adjusted_rand_index(truth[main_index], predicted[main_index]) else NA_real_,
    MainRetention = if (any(main_index)) mean(predicted[main_index] %in% seq_len(main_k)) else NA_real_,
    MinorityPrecision = precision,
    MinorityRecall = recall,
    MinorityF1 = f1,
    NoiseRecall = if (any(noise_index)) mean(predicted[noise_index] == 0L) else NA_real_,
    NoiseFalsePositiveRate = if (any(noise_index)) mean(predicted[noise_index] > main_k) else NA_real_,
    ExtraClusters = max(0L, length(unique(predicted[predicted > 0L])) - main_k),
    AnyMicroDiscovery = as.numeric(any(predicted > main_k)),
    row.names = NULL
  )
}
