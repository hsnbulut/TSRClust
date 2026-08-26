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
         eps_grid_size = 31L,
         standardize = TRUE, ...)
}
\arguments{
\item{x}{Numeric matrix or data frame.}
\item{k}{Number of dominant macro clusters.}
\item{alpha}{Trimming proportion used in the macro stage.}
\item{restr_fact}{TCLUST eigenvalue restriction factor for fixed macro fitting. When \code{auto_macro = TRUE}, omit this argument to scan \code{c_grid}; supply it only to constrain automatic selection to a single restriction factor.}
\item{eps}{Optional DBSCAN radius. If \code{NULL}, it is selected from the residual set.}
\item{min_pts}{Minimum number of residual points needed for a DBSCAN structure.}
\item{macro_method}{Macro-stage method: \code{"tclust"} or \code{"tkmeans"}.}
\item{auto_macro}{If \code{TRUE}, select TCLUST macro parameters over candidate grids. By default this scans \code{k_grid}, \code{alpha_grid}, and \code{c_grid}.}
\item{k_grid, alpha_grid, c_grid}{Candidate grids used when \code{auto_macro = TRUE}.}
\item{macro_ic}{Information criterion used by \code{tclustIC}.}
\item{validated}{If \code{TRUE}, validate residual structures against diffuse-noise reference samples.}
\item{null_reps}{Number of diffuse-noise reference samples.}
\item{structure_alpha}{Validation significance level for residual structure.}
\item{min_density_ratio}{Minimum density contrast required for residual structures.}
\item{novelty_prob}{Quantile used for macro-cluster novelty screening.}
\item{eps_quantile}{Nearest-neighbor distance quantile used when \code{validated = FALSE}.}
\item{eps_grid_size}{Number of quantile values in the residual epsilon candidate grid used when \code{validated = TRUE}. The manuscript experiments used 31.}
\item{standardize}{If \code{TRUE}, robustly standardize the analysis space before the macro stage; Stage-I clustering, residual discovery, and macro-reference construction then use the same standardized coordinates.}
\item{...}{Additional arguments passed to the macro-stage method.}
}
\value{
An object of class \code{tsrclust} containing final labels, raw macro labels,
residual indices, macro parameters, and micro-stage diagnostics.
}
\description{
Fits TSRClust by estimating dominant clusters in a trimmed robust macro-stage
and then applying density-based residual structural discovery only to the
trimmed observations. By default, the input is robustly standardized before
the macro stage so that the macro fit and residual density analysis are
performed in the same analysis space.
}
\details{
The package default \code{null_reps = 19L} is intended for interactive use.
The manuscript simulations and real-data analyses used \code{null_reps = 49L},
\code{structure_alpha = 0.05}, \code{min_density_ratio = 2},
\code{novelty_prob = 0.995}, \code{eps_quantile = 0.90}, and
\code{eps_grid_size = 31L}. These are tuning settings; the algorithmic
ordering is the same.

For \code{auto_macro = TRUE}, \code{tsrclust()} scans the supplied
\code{k_grid}, \code{alpha_grid}, and \code{c_grid} unless \code{k},
\code{alpha}, or \code{restr_fact} are explicitly supplied to constrain one of
those dimensions. The ordinary fixed-macro default remains
\code{auto_macro = FALSE, restr_fact = 20}.
}
\examples{
set.seed(1)
dat <- simulate_tsr_data(n = 200, scenario = "M1", seed = 1)
fit <- tsrclust(dat$x, k = 2, alpha = 0.10, validated = FALSE,
                nstart = 10)
fit
summary(fit)
}
