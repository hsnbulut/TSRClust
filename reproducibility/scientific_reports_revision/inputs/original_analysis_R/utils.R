`%||%` <- function(x, y) if (is.null(x) || length(x) == 0L || is.na(x[1L])) y else x

as_int <- function(x, default = NA_integer_) {
  if (length(x) == 0L || is.na(x) || x == "") return(default)
  as.integer(x)
}

as_num <- function(x, default = NA_real_) {
  if (length(x) == 0L || is.na(x) || x == "") return(default)
  as.numeric(x)
}

as_flag <- function(x) {
  if (is.logical(x)) return(isTRUE(x))
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

split_ints <- function(x) {
  if (length(x) == 0L || is.na(x) || x == "") return(integer(0L))
  as.integer(strsplit(as.character(x), ";", fixed = TRUE)[[1L]])
}

safe_dir_create <- function(path) {
  if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

safe_ratio <- function(num, den, undefined = NA_real_) {
  if (!is.finite(den) || den == 0) undefined else num / den
}

choose2 <- function(z) z * (z - 1) / 2

adjusted_rand_index <- function(truth, predicted) {
  truth <- as.vector(truth)
  predicted <- as.vector(predicted)
  if (length(truth) != length(predicted)) stop("ARI inputs have different lengths.")
  contingency <- table(truth, predicted)
  n <- sum(contingency)
  if (n < 2L) return(1)
  index <- sum(choose2(contingency))
  row_pairs <- sum(choose2(rowSums(contingency)))
  col_pairs <- sum(choose2(colSums(contingency)))
  total_pairs <- choose2(n)
  expected <- row_pairs * col_pairs / total_pairs
  maximum <- 0.5 * (row_pairs + col_pairs)
  denominator <- maximum - expected
  if (abs(denominator) <= sqrt(.Machine$double.eps)) {
    return(if (index == maximum) 1 else 0)
  }
  (index - expected) / denominator
}

robust_standardize <- function(x) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  center <- apply(x, 2L, stats::median)
  scale_value <- apply(x, 2L, stats::mad)
  fallback <- apply(x, 2L, stats::sd)
  bad <- !is.finite(scale_value) | scale_value <= sqrt(.Machine$double.eps)
  scale_value[bad] <- fallback[bad]
  scale_value[!is.finite(scale_value) | scale_value <= sqrt(.Machine$double.eps)] <- 1
  sweep(sweep(x, 2L, center, FUN = "-"), 2L, scale_value, FUN = "/")
}

seed_for <- function(base_seed, study_code, cell_id, replication, offset = 0L) {
  code <- sum(utf8ToInt(as.character(study_code))) %% 1000L
  value <- as.double(base_seed) + code * 1000000 + as.double(cell_id) * 10000 + replication + offset
  as.integer(value %% .Machine$integer.max)
}

run_safe <- function(expr) {
  warnings <- character(0L)
  started <- proc.time()[[3L]]
  value <- tryCatch(
    withCallingHandlers(
      expr,
      warning = function(w) {
        warnings <<- c(warnings, conditionMessage(w))
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) structure(list(message = conditionMessage(e)), class = "tsr_error")
  )
  elapsed <- proc.time()[[3L]] - started
  if (inherits(value, "tsr_error")) {
    list(ok = FALSE, value = NULL, error = value$message,
         warning = paste(unique(warnings), collapse = " | "), elapsed = elapsed)
  } else {
    list(ok = TRUE, value = value, error = "",
         warning = paste(unique(warnings), collapse = " | "), elapsed = elapsed)
  }
}

normalise_tkmeans_labels <- function(labels, k) {
  labels <- as.integer(labels)
  labels[labels == (as.integer(k) + 1L)] <- 0L
  labels
}

read_manifest_row <- function(manifest, task_id) {
  tab <- utils::read.csv(manifest, stringsAsFactors = FALSE, check.names = FALSE)
  row <- tab[tab$task_id == task_id, , drop = FALSE]
  if (nrow(row) != 1L) stop("Manifest task_id not found or duplicated: ", task_id)
  row
}
