# TSRClust

TSRClust implements Two-Stage Robust Clustering for Structural Minority
Discovery.

Stage I fits a trimmed robust macro-clustering model to protect the dominant
partition. When `standardize = TRUE`, robust standardisation is performed
before Stage-I macro clustering and the same standardised analysis space is
used for macro clustering, residual candidate discovery, density calculations
and macro-novelty assessment.

Stage II revisits only the trimmed residual set `R_0`. DBSCAN candidate search,
density scoring, macro-novelty checks, diffuse-noise null generation, and final
candidate validation all use the same analysis space used by Stage I.

## Installation

```r
install.packages("remotes")
remotes::install_github("hsnbulut/TSRClust")
```

From this local source checkout:

```r
install.packages("RPackageNew/TSRClust", repos = NULL, type = "source")
```

## Basic Example

```r
library(TSRClust)

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

summary(fit)
table(fit$cluster)
evaluate_tsrclust(dat$truth, fit$cluster)
```

The package default `null_reps = 19` is lighter for interactive use. The
manuscript simulations and real-data examples used `null_reps = 49`,
`structure_alpha = 0.05`, `min_density_ratio = 2`, `novelty_prob = 0.995`,
`eps_quantile = 0.90`, and `eps_grid_size = 31`.

When `auto_macro = TRUE`, `tsrclust()` scans the supplied `k_grid`,
`alpha_grid`, and `c_grid`. Supplying `k`, `alpha`, or `restr_fact` explicitly
constrains that dimension; otherwise the corresponding grid is used. The
ordinary fixed-macro default remains `auto_macro = FALSE` with
`restr_fact = 20`.

## Reproducibility

Final manuscript result files are organised under `results/`:

- `raw_results/`: replication-level simulation exports
- `summary_results/`: manuscript summary tables
- `real_data_results/`: olive-oil and Shuttle analysis outputs
- `simulation_design/`: final simulation manifests

The `reproducibility/` directory contains the analysis scripts used to generate
the simulations, real-data analyses, summaries, diagnostics, and manuscript
tables and figures. Temporary HPC worker files and obsolete exploratory scripts
are not included.

The package is not available on CRAN.

## Shuttle Data

The package does not redistribute the Statlog Shuttle data. The helper below
loads the data from the `mlbench` package when it is installed:

```r
shuttle <- load_shuttle_data()
str(shuttle)
```

The official UCI dataset DOI is `10.24432/C5WS31`.
