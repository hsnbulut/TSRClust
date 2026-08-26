\name{TSRClust-package}
\alias{TSRClust-package}
\alias{TSRClust}
\docType{package}
\title{Two-Stage Robust Clustering for Structural Minority Discovery}
\description{
The TSRClust package implements a two-stage robust clustering workflow. The
first stage estimates dominant clusters by a trimmed robust clustering method.
The second stage revisits the trimmed residual set and searches for dense
structural minorities while retaining diffuse residual observations as noise.
When \code{standardize = TRUE}, robust standardisation is performed before
Stage-I macro clustering, and the same standardised analysis space is used
for macro clustering, residual discovery, and macro novelty checks.
}
\author{
Hasan Bulut
}
