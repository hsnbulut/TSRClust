print.tsrclust <- function(x, ...) {
  cat("TSRClust fit\n")
  cat("  Macro method: ", x$macro_method, "\n", sep = "")
  cat("  Main clusters: ", x$n_main, "\n", sep = "")
  cat("  Structural minorities: ", x$n_minority, "\n", sep = "")
  cat("  Residual noise points: ", x$n_noise, "\n", sep = "")
  cat("  Standardized analysis space: ", x$standardized, "\n", sep = "")
  invisible(x)
}

summary.tsrclust <- function(object, ...) {
  labels <- object$cluster
  type <- ifelse(labels == 0L, "noise",
                 ifelse(labels <= object$n_main, "main", "structural_minority"))
  out <- list(
    call = object$call,
    macro_parameters = object$macro_parameters,
    micro_info = object$micro_info,
    counts = data.frame(
      label = as.integer(names(table(labels))),
      count = as.integer(table(labels)),
      type = tapply(type, labels, function(z) z[1L]),
      row.names = NULL
    )
  )
  class(out) <- "summary.tsrclust"
  out
}

print.summary.tsrclust <- function(x, ...) {
  cat("TSRClust summary\n")
  cat("Macro parameters:\n")
  print(x$macro_parameters)
  cat("\nCluster counts:\n")
  print(x$counts, row.names = FALSE)
  invisible(x)
}

plot.tsrclust <- function(x, data, pca = TRUE, ...) {
  data <- .as_numeric_matrix(data)
  if (nrow(data) != length(x$cluster)) {
    stop("data must have the same number of rows as the fitted object.", call. = FALSE)
  }
  if (ncol(data) > 2L && pca) {
    coords <- stats::prcomp(data, scale. = TRUE)$x[, 1L:2L, drop = FALSE]
    xlab <- "PC1"
    ylab <- "PC2"
  } else {
    if (ncol(data) < 2L) {
      stop("plot.tsrclust requires at least two variables.", call. = FALSE)
    }
    coords <- data[, 1L:2L, drop = FALSE]
    xlab <- colnames(data)[1L]
    ylab <- colnames(data)[2L]
  }
  type <- ifelse(x$cluster == 0L, "Noise",
                 ifelse(x$cluster <= x$n_main, "Main clusters", "Structural minorities"))

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    plot_df <- data.frame(X1 = coords[, 1L], X2 = coords[, 2L],
                          Cluster = factor(x$cluster),
                          Type = factor(type, levels = c("Main clusters", "Structural minorities", "Noise")))
    return(
      ggplot2::ggplot(plot_df, ggplot2::aes_string("X1", "X2", color = "Cluster", shape = "Type")) +
        ggplot2::geom_point(alpha = 0.75, size = 2) +
        ggplot2::labs(x = xlab, y = ylab, color = "Cluster", shape = "Type") +
        ggplot2::theme_minimal()
    )
  }

  graphics::plot(coords[, 1L], coords[, 2L], col = as.integer(factor(x$cluster)) + 1L,
                 pch = ifelse(type == "Noise", 4L, ifelse(type == "Structural minorities", 17L, 16L)),
                 xlab = xlab, ylab = ylab, ...)
  graphics::legend("topright", legend = unique(type), pch = c(16L, 17L, 4L), bty = "n")
  invisible(x)
}
