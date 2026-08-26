TSR_CONFIG <- list(
  base_seed = 20260805L,
  nstart_main = 500L,
  nstart_reference = 500L,
  null_reps_main = 49L,
  structure_alpha = 0.05,
  min_density_ratio = 2,
  novelty_prob = 0.995,
  eps_quantile = 0.90,
  eps_grid_size = 31L,
  block_size = 10L,
  tclust_niter1 = 3L,
  tclust_niter2 = 20L,
  tclust_nkeep = 5L,
  auto_k_grid = 1:5,
  auto_alpha_grid = c(0.05, 0.10, 0.15, 0.20),
  auto_c_grid = c(1, 5, 20, 100),
  auto_ic = "MIXMIX"
)
