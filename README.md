# TSRClust

`TSRClust` implements Two-Stage Robust Clustering for structural minority
discovery. The method first estimates dominant clusters using a trimmed robust
macro-stage, then revisits only the trimmed residual set to detect dense
minority structures while leaving diffuse residual observations as noise.

## Installation

After public release on GitHub:

```r
install.packages("remotes")
remotes::install_github("hsnbulut/TSRClust")
```

From a local source checkout:

```r
install.packages("TSRClust", repos = NULL, type = "source")
```

or, during development:

```r
pkgload::load_all("TSRClust")
```

## Basic Example

```r
library(TSRClust)

set.seed(1)
dat <- simulate_tsr_data(n = 500, scenario = "M1")

fit <- tsrclust(
  dat$x,
  k = 2,
  alpha = 0.10,
  restr_fact = 20
)

fit
summary(fit)
evaluate_tsrclust(dat$truth, fit$cluster)
plot(fit, dat$x)
```

## Shuttle Data Example

The package does not redistribute the Statlog Shuttle data. The helper below
loads the data from the `mlbench` package when it is installed:

```r
shuttle <- load_shuttle_data()
str(shuttle)
```

This keeps the package lightweight while still providing a real-data example
for demonstrations and teaching.
