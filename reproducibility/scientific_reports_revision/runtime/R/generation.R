revision_data <- function(row,seed) {
  # Reuse submitted generators; only additional background mechanisms are new.
  z <- row
  extra <- row$noise_type %in% c('correlated_gaussian','heterogeneous','heavy_tail')
  if(extra) z$noise_type <- 'uniform'
  d <- legacy$simulate_protocol_data(z,seed)
  if(extra) {
    set.seed(seed+10000000L) # reserved data-generation substream, no algorithm/null overlap
    idx <- which(d$truth==0L); n <- length(idx); p <- row$p
    if(row$noise_type=='correlated_gaussian') {
      S <- matrix(.8,p,p);diag(S)<-1
      d$x[idx,] <- legacy$rmvnorm_chol(n,rep(0,p),S*(.75*row$separation)^2)
    }
    if(row$noise_type=='heterogeneous') d$x[idx,] <- matrix((rbeta(n*p,2,1)*2-1)*.8*row$separation,n,p)
    if(row$noise_type=='heavy_tail') d$x[idx,] <- legacy$rmvt_chol(n,rep(0,p),diag((.75*row$separation)^2,p),3)
  }
  if(row$mode=='residual_only') {
    set.seed(seed)
    d$x <- matrix(runif(row$nR*row$p,-4,4),row$nR,row$p)
    d$truth <- integer(row$nR); d$main_k <- 0L; d$q <- 0L
  }
  d
}
