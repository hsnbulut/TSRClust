\name{simulate_tsr_data}
\alias{simulate_tsr_data}
\title{Simulate TSRClust Benchmark Data}
\usage{
simulate_tsr_data(n = 500, p = 2,
                  scenario = c("N0", "M1", "M2", "M3", "M4"),
                  structure = c("spherical", "elliptical"),
                  separation = 6, seed = NULL)
}
\arguments{
\item{n}{Sample size.}
\item{p}{Number of variables.}
\item{scenario}{Residual-pool scenario.}
\item{structure}{Dominant-cluster covariance structure.}
\item{separation}{Distance between dominant cluster centers.}
\item{seed}{Optional random seed.}
}
\value{
A list with \code{x}, \code{truth}, scenario metadata, and the true number of
dominant and minority clusters.
}
\description{
Generates the simulation scenarios used in the TSRClust manuscript. The data
contain two dominant Gaussian clusters and a residual pool consisting of pure
noise, structural minorities, or both.
}
\examples{
dat <- simulate_tsr_data(n = 100, scenario = "M2", seed = 1)
table(dat$truth)
}
