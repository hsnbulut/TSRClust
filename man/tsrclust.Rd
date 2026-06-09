\name{tsrclust}
\alias{tsrclust}
\title{Two-Stage Robust Clustering}
\usage{
tsrclust(x, k = NULL, alpha = NULL, restr_fact = 20,
         eps = NULL, min_pts = NULL,
         macro_method = c("tclust", "tkmeans"),
         auto_macro = FALSE, k_grid = 1:5,
         alpha_grid = seq(0, 0.20, length.out = 6),
         c_grid = c(1, 2, 4, 8, 16, 32, 64, 128),
         macro_ic = c("MIXMIX", "MIXCLA", "CLACLA"),
         validated = TRUE, null_reps = 19L,
         structure_alpha = 0.05, min_density_ratio = 2,
         novelty_prob = 0.995, eps_quantile = 0.90,
         standardize = TRUE, ...)
}
\arguments{
\item{x}{Numeric matrix or data frame.}
\item{k}{Number of dominant macro clusters.}
\item{alpha}{Trimming proportion used in the macro stage.}
\item{restr_fact}{TCLUST eigenvalue restriction factor.}
\item{eps}{Optional DBSCAN radius. If \code{NULL}, it is selected from the residual set.}
\item{min_pts}{Minimum number of residual points needed for a DBSCAN structure.}
\item{macro_method}{Macro-stage method: \code{"tclust"} or \code{"tkmeans"}.}
\item{auto_macro}{If \code{TRUE}, select TCLUST macro parameters over candidate grids.}
\item{k_grid, alpha_grid, c_grid}{Candidate grids used when \code{auto_macro = TRUE}.}
\item{macro_ic}{Information criterion used by \code{tclustIC}.}
\item{validated}{If \code{TRUE}, validate residual structures against diffuse-noise reference samples.}
\item{null_reps}{Number of diffuse-noise reference samples.}
\item{structure_alpha}{Validation significance level for residual structure.}
\item{min_density_ratio}{Minimum density contrast required for residual structures.}
\item{novelty_prob}{Quantile used for macro-cluster novelty screening.}
\item{eps_quantile}{Nearest-neighbor distance quantile used when \code{validated = FALSE}.}
\item{standardize}{If \code{TRUE}, robustly standardize the DBSCAN distance space.}
\item{...}{Additional arguments passed to the macro-stage method.}
}
\value{
An object of class \code{tsrclust} containing final labels, raw macro labels,
residual indices, macro parameters, and micro-stage diagnostics.
}
\description{
Fits TSRClust by estimating dominant clusters in a trimmed robust macro-stage
and then applying density-based residual structural discovery only to the
trimmed observations.
}
\examples{
set.seed(1)
dat <- simulate_tsr_data(n = 200, scenario = "M1")
fit <- tsrclust(dat$x, k = 2, alpha = 0.10, validated = FALSE)
fit
summary(fit)
}
