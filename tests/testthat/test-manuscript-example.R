test_that("manuscript example executes", {
  set.seed(1)
  dat <- simulate_tsr_data(n = 500, p = 2, scenario = "M1",
                           structure = "elliptical", seed = 1)
  fit <- tsrclust(dat$x, k = 2, alpha = 0.10, restr_fact = 20,
                  validated = TRUE, null_reps = 49,
                  structure_alpha = 0.05,
                  min_density_ratio = 2,
                  novelty_prob = 0.995,
                  eps_grid_size = 31,
                  nstart = 500)

  expect_s3_class(fit, "tsrclust")
  expect_equal(length(fit$cluster), nrow(dat$x))
  expect_s3_class(summary(fit), "summary.tsrclust")
  expect_s3_class(evaluate_tsrclust(dat$truth, fit$cluster), "data.frame")
})
