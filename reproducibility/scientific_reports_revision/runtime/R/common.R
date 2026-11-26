bind_rows <- function(xs) {
  if (!length(xs)) return(data.frame())
  nms <- unique(unlist(lapply(xs,names)))
  do.call(rbind,lapply(xs,function(x) {for(nm in setdiff(nms,names(x))) x[[nm]] <- NA; x[nms]}))
}
seed_id <- function(dataset, rep, stream, repeat_id=1L) {
  stopifnot(dataset>0,dataset<1000,rep>0,rep<10000,stream %in% 0:2,repeat_id %in% 1:9)
  # Disjoint namespaces: 0/1/2=data/algorithm/null, disjoint repeats.
  as.integer(stream*600000000 + (repeat_id-1)*60000000 + dataset*10000 + rep)
}
atomic_rds <- function(object,path) {
  dir.create(dirname(path),recursive=TRUE,showWarnings=FALSE)
  tmp <- paste0(path,'.',Sys.getpid(),'.tmp')
  on.exit(unlink(tmp),add=TRUE)
  saveRDS(object,tmp,compress='gzip')
  if (!file.rename(tmp,path)) stop('atomic rename failed: ',path)
}
exact_interval <- function(k,n) {
  if (!n) return(c(NA_real_,NA_real_))
  unname(binom.test(k,n)$conf.int)
}
check_labels <- function(y,pred) {
  if(length(y)!=length(pred)||anyNA(pred)||any(!is.finite(pred))||any(pred<0)||any(pred!=as.integer(pred))) stop('Invalid partition')
}
