\name{load_shuttle_data}
\alias{load_shuttle_data}
\title{Load the Statlog Shuttle Data}
\usage{
load_shuttle_data(scale = TRUE, drop_time = FALSE)
}
\arguments{
\item{scale}{If \code{TRUE}, standardize the numeric variables.}
\item{drop_time}{If \code{TRUE}, remove the first numeric variable, which is time.}
}
\value{
A list with numeric matrix \code{x}, factor vector \code{truth}, and source
metadata.
}
\description{
Loads the Statlog Shuttle data from the \code{mlbench} package and returns it
in a form suitable for TSRClust examples. The data are not redistributed by
TSRClust; users should install \code{mlbench} to access the original dataset.
}
\references{
UCI Machine Learning Repository, Statlog Shuttle data set,
\doi{10.24432/C5WS31}.
}
\examples{
if (requireNamespace("mlbench", quietly = TRUE)) {
  shuttle <- load_shuttle_data()
  table(shuttle$truth)
}
}
