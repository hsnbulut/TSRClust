tsr_scenarios <- function() {
  data.frame(
    scenario = c("N0", "M1", "M2", "M3", "M4"),
    description = c(
      "pure noise",
      "one compact minority plus noise",
      "one smaller compact minority plus noise",
      "one diffuse minority plus noise",
      "two compact minorities plus noise"
    ),
    minority_ratio = c(0.000, 0.050, 0.025, 0.050, 0.050),
    noise_ratio = c(0.100, 0.050, 0.075, 0.050, 0.050),
    n_minority_clusters = c(0L, 1L, 1L, 1L, 2L),
    minority_sd = c(NA, 0.20, 0.20, 0.55, 0.20),
    stringsAsFactors = FALSE
  )
}

simulate_tsr_data <- function(n = 500, p = 2,
                              scenario = c("N0", "M1", "M2", "M3", "M4"),
                              structure = c("spherical", "elliptical"),
                              separation = 6,
                              seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  scenario <- match.arg(scenario)
  structure <- match.arg(structure)
  catalog <- tsr_scenarios()
  scenario_row <- catalog[catalog$scenario == scenario, ]

  n_noise <- round(n * scenario_row$noise_ratio)
  n_minority <- round(n * scenario_row$minority_ratio)
  n_main <- n - n_noise - n_minority
  main_counts <- .allocate_counts(n_main, 2L)
  main_means <- list(rep(0, p), rep(0, p))
  main_means[[1L]][1L] <- -separation / 2
  main_means[[2L]][1L] <- separation / 2

  x_main <- do.call(rbind, lapply(seq_len(2L), function(g) {
    .rmvnorm(main_counts[g], main_means[[g]], .main_covariance(p, structure, g))
  }))
  y_main <- rep(seq_len(2L), main_counts)

  q <- scenario_row$n_minority_clusters
  if (q > 0L) {
    minority_counts <- .allocate_counts(n_minority, q)
    minority_centers <- lapply(seq_len(q), function(g) {
      center <- rep(0, p)
      center[2L] <- if (q == 1L) 0.75 * separation else c(0.75, -0.75)[g] * separation
      center
    })
    minority_sigma <- diag(scenario_row$minority_sd^2, p)
    x_minority <- do.call(rbind, lapply(seq_len(q), function(g) {
      .rmvnorm(minority_counts[g], minority_centers[[g]], minority_sigma)
    }))
    y_minority <- rep(2L + seq_len(q), minority_counts)
  } else {
    x_minority <- matrix(numeric(0L), nrow = 0L, ncol = p)
    y_minority <- integer(0L)
  }

  noise_limits <- matrix(rep(c(-separation, separation), p), ncol = 2L, byrow = TRUE)
  if (p > 2L) {
    noise_limits[3L:p, ] <- matrix(rep(c(-4, 4), p - 2L), ncol = 2L, byrow = TRUE)
  }
  x_noise <- vapply(seq_len(p), function(j) {
    stats::runif(n_noise, min = noise_limits[j, 1L], max = noise_limits[j, 2L])
  }, numeric(n_noise))
  x_noise <- matrix(x_noise, ncol = p)
  y_noise <- rep(0L, n_noise)

  x <- rbind(x_main, x_minority, x_noise)
  truth <- c(y_main, y_minority, y_noise)
  permutation <- sample.int(nrow(x))
  x <- x[permutation, , drop = FALSE]
  truth <- truth[permutation]
  colnames(x) <- paste0("V", seq_len(p))

  out <- list(
    x = x,
    truth = truth,
    scenario = scenario,
    main_k = 2L,
    n_minority_clusters = q,
    minority_ratio = scenario_row$minority_ratio,
    noise_ratio = scenario_row$noise_ratio
  )
  class(out) <- "tsr_data"
  out
}

.allocate_counts <- function(total, groups) {
  counts <- rep(total %/% groups, groups)
  if (total %% groups > 0L) {
    counts[seq_len(total %% groups)] <- counts[seq_len(total %% groups)] + 1L
  }
  counts
}

.main_covariance <- function(p, structure, cluster_id) {
  sigma <- diag(p)
  if (structure == "elliptical" && p >= 2L) {
    if (cluster_id == 1L) {
      sigma[1L:2L, 1L:2L] <- matrix(c(1.8, 0.9, 0.9, 1.0), nrow = 2L)
    } else {
      sigma[1L:2L, 1L:2L] <- matrix(c(1.0, -0.6, -0.6, 1.8), nrow = 2L)
    }
  }
  sigma
}

.rmvnorm <- function(n, mean, sigma) {
  if (n <= 0L) {
    return(matrix(numeric(0L), nrow = 0L, ncol = length(mean)))
  }
  z <- matrix(stats::rnorm(n * length(mean)), nrow = n)
  sweep(z %*% chol(sigma), 2L, mean, FUN = "+")
}
