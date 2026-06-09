load_shuttle_data <- function(scale = TRUE, drop_time = FALSE) {
  if (!requireNamespace("mlbench", quietly = TRUE)) {
    stop("Package 'mlbench' is required to load the Shuttle data.", call. = FALSE)
  }

  utils::data("Shuttle", package = "mlbench", envir = environment())
  shuttle <- get("Shuttle", envir = environment())

  class_column <- which(vapply(shuttle, is.factor, logical(1L)))
  if (length(class_column) != 1L) {
    stop("Could not uniquely identify the Shuttle class column.", call. = FALSE)
  }

  x <- shuttle[, -class_column, drop = FALSE]
  truth <- droplevels(shuttle[[class_column]])
  numeric_columns <- vapply(x, is.numeric, logical(1L))
  x <- x[, numeric_columns, drop = FALSE]

  if (drop_time && ncol(x) > 1L) {
    x <- x[, -1L, drop = FALSE]
  }

  non_constant <- vapply(x, function(z) stats::sd(z) > 0, logical(1L))
  x <- x[, non_constant, drop = FALSE]
  x <- as.matrix(x)
  if (scale) {
    x <- scale(x)
  }

  list(
    x = x,
    truth = truth,
    source = "Statlog Shuttle data from mlbench::Shuttle",
    citation = "UCI Machine Learning Repository, doi:10.24432/C5WS31"
  )
}
