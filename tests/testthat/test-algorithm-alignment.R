test_that("standardize = TRUE sends standardized data to Stage-I TCLUST", {
  dat <- simulate_tsr_data(n = 160, p = 2, scenario = "M1", seed = 11)
  x <- dat$x
  x[, 2] <- 100 * x[, 2] + 500

  set.seed(101)
  direct <- tclust::tclust(TSRClust:::.robust_standardize(x), k = 2,
                           alpha = 0.10, restr.fact = 20, nstart = 20)
  set.seed(101)
  fit <- tsrclust(x, k = 2, alpha = 0.10, restr_fact = 20,
                  standardize = TRUE, validated = FALSE, nstart = 20)

  expect_equal(fit$raw_cluster, as.integer(direct$cluster))
  expect_true(fit$standardized)
})

test_that("standardize = FALSE leaves Stage-I analysis space unchanged", {
  dat <- simulate_tsr_data(n = 160, p = 2, scenario = "M1", seed = 12)
  x <- dat$x
  x[, 2] <- 100 * x[, 2] + 500

  set.seed(102)
  direct <- tclust::tclust(x, k = 2, alpha = 0.10,
                           restr.fact = 20, nstart = 20)
  set.seed(102)
  fit <- tsrclust(x, k = 2, alpha = 0.10, restr_fact = 20,
                  standardize = FALSE, validated = FALSE, nstart = 20)

  expect_equal(fit$raw_cluster, as.integer(direct$cluster))
  expect_false(fit$standardized)
})

test_that("Stage II uses the same standardized analysis matrix", {
  dat <- simulate_tsr_data(n = 180, p = 2, scenario = "M1", seed = 13)
  x <- dat$x
  x[, 2] <- 25 * x[, 2] - 100

  set.seed(103)
  fit <- tsrclust(x, k = 2, alpha = 0.10, restr_fact = 20,
                  standardize = TRUE, validated = FALSE, nstart = 20)

  x_analysis <- TSRClust:::.robust_standardize(x)
  residual_x <- x_analysis[fit$residual_index, , drop = FALSE]
  eps <- TSRClust:::.select_eps_knn(residual_x, fit$min_pts, 0.90)
  direct_micro <- if (is.na(eps)) {
    integer(nrow(residual_x))
  } else {
    as.integer(dbscan::dbscan(residual_x, eps = eps, minPts = fit$min_pts)$cluster)
  }

  expect_equal(fit$micro_cluster, direct_micro)
})

test_that("known simulated M1 data produce a valid tsrclust object", {
  dat <- simulate_tsr_data(n = 200, p = 2, scenario = "M1", seed = 14)
  set.seed(104)
  fit <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                  validated = FALSE, nstart = 20)

  expect_s3_class(fit, "tsrclust")
  expect_equal(length(fit$cluster), nrow(dat$x))
  expect_true(all(fit$cluster >= 0))
})

test_that("cluster labels and residual handling are valid", {
  dat <- simulate_tsr_data(n = 180, p = 2, scenario = "N0", seed = 15)
  set.seed(105)
  fit <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                  validated = FALSE, nstart = 20)

  expect_equal(which(fit$raw_cluster == 0L), fit$residual_index)
  expect_equal(length(fit$micro_cluster), length(fit$residual_index))
  expect_true(all(fit$cluster[fit$cluster > fit$n_main] > fit$n_main))
})

test_that("validated = FALSE and validated = TRUE both work", {
  dat <- simulate_tsr_data(n = 180, p = 2, scenario = "M1", seed = 16)

  set.seed(106)
  unvalidated <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                          validated = FALSE, nstart = 20)
  set.seed(106)
  validated <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                        validated = TRUE, null_reps = 3,
                        eps_grid_size = 5, nstart = 20)

  expect_s3_class(unvalidated, "tsrclust")
  expect_s3_class(validated, "tsrclust")
  expect_true(is.list(validated$micro_info))
})

test_that("automatic macro selection still works after preprocessing", {
  dat <- simulate_tsr_data(n = 140, p = 2, scenario = "M1", seed = 17)

  set.seed(107)
  fit <- tsrclust(dat$x, auto_macro = TRUE,
                  k_grid = 2, alpha_grid = 0.10, c_grid = 20,
                  validated = FALSE, nstart = 10)

  expect_s3_class(fit, "tsrclust")
  expect_equal(fit$n_main, 2L)
  expect_equal(fit$macro_parameters$alpha, 0.10)
})

test_that("automatic macro selection scans c_grid when restr_fact is omitted", {
  dat <- simulate_tsr_data(n = 140, p = 2, scenario = "M1", seed = 171)

  set.seed(117)
  fit <- tsrclust(dat$x, auto_macro = TRUE,
                  k_grid = 2, alpha_grid = 0.10, c_grid = c(1, 20),
                  validated = FALSE, nstart = 10)

  expect_s3_class(fit, "tsrclust")
  expect_equal(fit$macro_selection$scanned_c_grid, c(1, 20))

  set.seed(117)
  constrained <- tsrclust(dat$x, auto_macro = TRUE, restr_fact = 20,
                          k_grid = 2, alpha_grid = 0.10, c_grid = c(1, 20),
                          validated = FALSE, nstart = 10)

  expect_equal(constrained$macro_selection$scanned_c_grid, 20)
})

test_that("print, summary, plot, and evaluate_tsrclust do not break", {
  dat <- simulate_tsr_data(n = 160, p = 2, scenario = "M1", seed = 18)
  set.seed(108)
  fit <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                  validated = FALSE, nstart = 20)

  expect_output(print(fit), "TSRClust fit")
  expect_s3_class(summary(fit), "summary.tsrclust")
  expect_silent(plot(fit, dat$x))
  expect_s3_class(evaluate_tsrclust(dat$truth, fit$cluster), "data.frame")
})

test_that("fixed-seed runs are reproducible", {
  dat <- simulate_tsr_data(n = 160, p = 2, scenario = "M1", seed = 19)

  set.seed(109)
  fit1 <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                   validated = TRUE, null_reps = 3,
                   eps_grid_size = 5, nstart = 20)
  set.seed(109)
  fit2 <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                   validated = TRUE, null_reps = 3,
                   eps_grid_size = 5, nstart = 20)

  expect_equal(fit1$raw_cluster, fit2$raw_cluster)
  expect_equal(fit1$cluster, fit2$cluster)
  expect_equal(fit1$micro_info$p_value, fit2$micro_info$p_value)
})
