\name{evaluate_tsrclust}
\alias{evaluate_tsrclust}
\title{Evaluate a TSRClust Partition}
\usage{
evaluate_tsrclust(truth, predicted, main_k = 2,
                  predicted_positive = predicted > main_k)
}
\arguments{
\item{truth}{Integer vector of true labels. Label 0 denotes noise.}
\item{predicted}{Integer vector of predicted labels. Label 0 denotes residual noise.}
\item{main_k}{Number of dominant clusters.}
\item{predicted_positive}{Logical vector indicating observations predicted as structural minorities.}
}
\value{
A one-row data frame with global ARI, main-cluster ARI, minority precision,
minority recall, minority F1, noise recall, and related summary measures.
}
\description{
Computes external evaluation metrics used in the TSRClust simulation study.
}
\examples{
dat <- simulate_tsr_data(n = 100, scenario = "M1", seed = 1)
fit <- tsrclust(dat$x, k = 2, alpha = 0.10, validated = FALSE)
evaluate_tsrclust(dat$truth, fit$cluster)
}
