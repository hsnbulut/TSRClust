fit_otrimle_fixed <- function(x, G, erc = 20, npr.max = 0.50) {
  fit <- otrimle::otrimle(
    data = x,
    G = as.integer(G),
    npr.max = npr.max,
    erc = erc,
    ncores = 1,
    monitor = FALSE
  )
  if (is.null(fit$cluster)) stop("OTRIMLE returned no clustering.")
  list(cluster = as.integer(fit$cluster), code = fit$code,
       criterion = fit$criterion, noise_proportion = mean(fit$cluster == 0L), fit = fit)
}

fit_teigen_bic <- function(x, Gs) {
  fit <- teigen::teigen(
    x = x,
    Gs = as.integer(Gs),
    models = "all",
    init = "kmeans",
    scale = FALSE,
    verbose = FALSE,
    parallel.cores = FALSE
  )
  if (is.null(fit$classification)) stop("tEIGEN returned no classification.")
  list(cluster = as.integer(fit$classification), selected_G = fit$G,
       modelname = fit$modelname, bic = fit$bic, fit = fit)
}

# End-to-end runtime is defined as preprocessing + every stage needed to obtain final labels.
# For TSRClust, the macro fit is shared exactly with TCLUST; therefore total time is
# preprocessing + macro elapsed + incremental residual-stage elapsed.
