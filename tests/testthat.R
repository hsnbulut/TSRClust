if (requireNamespace("testthat", quietly = TRUE)) {
  library(testthat)
  library(TSRClust)
  test_check("TSRClust")
} else {
  message("testthat is not installed; skipping testthat suite.")
}
