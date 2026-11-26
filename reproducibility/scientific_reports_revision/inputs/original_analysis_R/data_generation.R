allocate_counts <- function(total, weights) {
  weights <- weights / sum(weights)
  raw <- total * weights
  counts <- floor(raw)
  remainder <- total - sum(counts)
  if (remainder > 0L) {
    order_idx <- order(raw - counts, decreasing = TRUE)
    counts[order_idx[seq_len(remainder)]] <- counts[order_idx[seq_len(remainder)]] + 1L
  }
  as.integer(counts)
}

regular_simplex_centers <- function(k, p, separation) {
  if (k == 1L) return(matrix(0, nrow = 1L, ncol = p))
  d <- k - 1L
  if (p < d) stop("p must be at least k-1 for the simplex dominant-cluster design.")
  H <- diag(k) - matrix(1 / k, k, k)
  eig <- eigen(H, symmetric = TRUE)
  coords <- eig$vectors[, seq_len(d), drop = FALSE] %*%
    diag(sqrt(pmax(eig$values[seq_len(d)], 0)), d, d)
  current <- sqrt(sum((coords[1L, ] - coords[2L, ])^2))
  coords <- coords * (separation / current)
  out <- matrix(0, nrow = k, ncol = p)
  out[, seq_len(d)] <- coords
  out
}

main_covariance <- function(p, structure, cluster_id) {
  if (structure == "spherical") return(diag(p))
  sigma <- diag(p)
  if (p >= 2L) {
    sign <- if (cluster_id %% 2L == 1L) 1 else -1
    sigma[1:2, 1:2] <- matrix(c(1.8, 0.75 * sign, 0.75 * sign, 0.8), 2L, 2L)
  }
  if (p >= 3L) sigma[3L, 3L] <- if (cluster_id %% 2L == 1L) 1.4 else 0.7
  sigma
}

rmvnorm_chol <- function(n, mean, sigma) {
  if (n <= 0L) return(matrix(numeric(0L), 0L, length(mean)))
  ridge <- max(mean(diag(sigma)), 1) * 1e-10
  L <- chol(sigma + diag(ridge, ncol(sigma)))
  z <- matrix(stats::rnorm(n * length(mean)), nrow = n)
  sweep(z %*% L, 2L, mean, FUN = "+")
}

rmvt_chol <- function(n, mean, sigma, df) {
  if (n <= 0L) return(matrix(numeric(0L), 0L, length(mean)))
  z <- rmvnorm_chol(n, rep(0, length(mean)), sigma)
  scale <- sqrt(stats::rchisq(n, df = df) / df)
  sweep(z, 1L, scale, FUN = "/") + matrix(mean, nrow = n, ncol = length(mean), byrow = TRUE)
}

main_weights_vector <- function(k, mode) {
  if (mode == "balanced") return(rep(1 / k, k))
  if (k == 2L) return(c(0.75, 0.25))
  if (k == 4L) return(c(0.55, 0.25, 0.15, 0.05))
  w <- rev(seq_len(k))
  w / sum(w)
}

minority_centers <- function(q, p, main_centers, separation, proximity) {
  if (q <= 0L) return(list())
  base <- colMeans(main_centers)
  radius <- switch(proximity,
                   isolated = 0.75 * separation,
                   adjacent = 2.5,
                   overlapping = 1.25,
                   0.75 * separation)
  centers <- vector("list", q)
  if (proximity %in% c("adjacent", "overlapping")) base <- main_centers[1L, ]
  for (g in seq_len(q)) {
    angle <- 2 * pi * (g - 1L) / max(q, 1L)
    c0 <- base
    if (p >= 2L) {
      c0[1L] <- c0[1L] + radius * cos(angle + pi / 2)
      c0[2L] <- c0[2L] + radius * sin(angle + pi / 2)
    } else {
      c0[1L] <- c0[1L] + if (g %% 2L) radius else -radius
    }
    centers[[g]] <- c0
  }
  centers
}

generate_minority_group <- function(n, center, shape, sd_value) {
  p <- length(center)
  if (n <= 0L) return(matrix(numeric(0L), 0L, p))
  if (shape == "gaussian") return(rmvnorm_chol(n, center, diag(sd_value^2, p)))
  if (shape == "elongated") {
    sigma <- diag(sd_value^2, p)
    sigma[1L, 1L] <- (3 * sd_value)^2
    if (p >= 2L) sigma[2L, 2L] <- (0.5 * sd_value)^2
    return(rmvnorm_chol(n, center, sigma))
  }
  if (shape == "ring") {
    theta <- stats::runif(n, 0, 2 * pi)
    radius <- pmax(0.05, stats::rnorm(n, mean = 1.5, sd = max(0.08, sd_value / 3)))
    x <- matrix(0, nrow = n, ncol = p)
    x[, 1L] <- center[1L] + radius * cos(theta)
    if (p >= 2L) x[, 2L] <- center[2L] + radius * sin(theta)
    if (p >= 3L) x[, 3L:p] <- matrix(stats::rnorm(n * (p - 2L), sd = sd_value), nrow = n)
    return(x)
  }
  stop("Unknown minority shape: ", shape)
}

generate_noise <- function(n, p, type, main_centers, minority_centers_list, separation) {
  if (n <= 0L) return(matrix(numeric(0L), 0L, p))
  all_centers <- main_centers
  if (length(minority_centers_list) > 0L) all_centers <- rbind(all_centers, do.call(rbind, minority_centers_list))
  lower <- apply(all_centers, 2L, min) - max(2.5, 0.35 * separation)
  upper <- apply(all_centers, 2L, max) + max(2.5, 0.35 * separation)
  if (type == "uniform") {
    x <- vapply(seq_len(p), function(j) stats::runif(n, lower[j], upper[j]), numeric(n))
    return(matrix(x, nrow = n, ncol = p))
  }
  if (type == "gaussian_diffuse") {
    return(rmvnorm_chol(n, colMeans(main_centers), diag((0.75 * separation)^2, p)))
  }
  if (type == "bridge") {
    a <- main_centers[1L, ]
    b <- main_centers[min(2L, nrow(main_centers)), ]
    t <- stats::runif(n)
    x <- (1 - t) * matrix(a, n, p, byrow = TRUE) + t * matrix(b, n, p, byrow = TRUE)
    x <- x + matrix(stats::rnorm(n * p, sd = 0.25), nrow = n)
    return(x)
  }
  if (type == "correlated_sheet") {
    direction <- rep(0, p)
    direction[1L] <- 1
    if (p >= 2L) direction[2L] <- 0.6
    direction <- direction / sqrt(sum(direction^2))
    z <- stats::runif(n, -0.8 * separation, 0.8 * separation)
    x <- outer(z, direction) + matrix(stats::rnorm(n * p, sd = 0.18), nrow = n)
    x <- sweep(x, 2L, colMeans(main_centers), FUN = "+")
    return(x)
  }
  stop("Unknown noise type: ", type)
}

simulate_protocol_data <- function(row, seed) {
  set.seed(seed)
  n <- as_int(row$n)
  p <- as_int(row$p)
  k0 <- as_int(row$main_k)
  separation <- as_num(row$separation)
  residual_total <- as_num(row$residual_total)
  minority_ratio <- as_num(row$minority_ratio)
  q <- as_int(row$q, 0L)
  noise_ratio <- max(0, residual_total - minority_ratio)
  main_total <- n - round(n * residual_total)
  minority_total <- round(n * minority_ratio)
  noise_total <- n - main_total - minority_total
  centers <- regular_simplex_centers(k0, p, separation)
  weights <- main_weights_vector(k0, as.character(row$main_weights))
  counts <- allocate_counts(main_total, weights)
  main_x <- vector("list", k0)
  for (g in seq_len(k0)) {
    sigma <- main_covariance(p, as.character(row$main_structure), g)
    dist <- as.character(row$main_distribution)
    if (dist == "normal") main_x[[g]] <- rmvnorm_chol(counts[g], centers[g, ], sigma)
    if (dist == "t5") main_x[[g]] <- rmvt_chol(counts[g], centers[g, ], sigma, 5)
    if (dist == "t3") main_x[[g]] <- rmvt_chol(counts[g], centers[g, ], sigma, 3)
  }
  x_main <- do.call(rbind, main_x)
  y_main <- rep(seq_len(k0), counts)
  m_centers <- minority_centers(q, p, centers, separation, as.character(row$minority_proximity))
  if (q > 0L && minority_total > 0L) {
    m_counts <- allocate_counts(minority_total, rep(1 / q, q))
    m_x <- lapply(seq_len(q), function(g) {
      generate_minority_group(m_counts[g], m_centers[[g]],
                              as.character(row$minority_shape), as_num(row$minority_sd))
    })
    x_min <- do.call(rbind, m_x)
    y_min <- rep(k0 + seq_len(q), m_counts)
  } else {
    x_min <- matrix(numeric(0L), 0L, p)
    y_min <- integer(0L)
  }
  x_noise <- generate_noise(noise_total, p, as.character(row$noise_type),
                            centers, m_centers, separation)
  y_noise <- rep(0L, noise_total)
  x <- rbind(x_main, x_min, x_noise)
  truth <- c(y_main, y_min, y_noise)
  perm <- sample.int(nrow(x))
  x <- x[perm, , drop = FALSE]
  truth <- truth[perm]
  colnames(x) <- paste0("V", seq_len(p))
  list(x = x, truth = truth, main_k = k0, q = q,
       noise_ratio = noise_ratio, minority_ratio = minority_ratio,
       residual_total = residual_total)
}
