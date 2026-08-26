source("config/config.R")
source("R/utils.R")
source("R/core_tsr.R")
source("R/evaluation_metrics.R")
source("R/competitors.R")
source("v2/config_v2.R")
source("v2/R/competitors_v2.R")

dir.create("v2/results/realdata", recursive=TRUE, showWarnings=FALSE)

simple_ari <- function(truth, pred) adjusted_rand_index(as.integer(factor(truth)), as.integer(factor(pred)))
cluster_purity <- function(truth, pred, positive_only=FALSE) {
  ids <- sort(unique(pred[if (positive_only) pred > 0 else rep(TRUE,length(pred))]))
  if (!length(ids)) return(data.frame())
  do.call(rbind,lapply(ids,function(g){
    idx <- which(pred==g); tt <- table(truth[idx]); best <- names(tt)[which.max(tt)]
    data.frame(cluster=g,size=length(idx),dominant_class=best,purity=max(tt)/length(idx),stringsAsFactors=FALSE)
  }))
}

# Orthonormal ILR coordinates without an external compositional-data dependency.
ilr_transform <- function(raw, pseudocount=1) {
  z <- as.matrix(raw) + pseudocount
  comp <- z / rowSums(z)
  D <- ncol(comp)
  A <- rbind(diag(D-1L), rep(-1, D-1L))
  B <- qr.Q(qr(A))
  log(comp) %*% B
}

choose_c_fixed_k <- function(x, k, alpha, cgrid) {
  ic <- tclust::tclustIC(x, kk=as.integer(k), cc=cgrid, alpha=alpha, whichIC="MIXMIX",
                         nstart=TSR_V2$nstart, niter1=TSR_CONFIG$tclust_niter1,
                         niter2=TSR_CONFIG$tclust_niter2, nkeep=TSR_CONFIG$tclust_nkeep,
                         parallel=FALSE, trace=FALSE, store_x=FALSE)
  v <- ic$MIXMIX; v[!is.finite(v)] <- Inf
  pos <- which(v==min(v),arr.ind=TRUE)[1,]
  as.numeric(ic$cc[pos[2]])
}

analyse_one <- function(name, x_raw, macro_truth, local_truth, k, alpha, c,
                        run_global_robust=TRUE) {
  tprep <- proc.time()[[3L]]
  x <- robust_standardize(x_raw)
  prep <- proc.time()[[3L]] - tprep

  set.seed(TSR_V2$base_seed + 100L)
  macro <- run_safe(fit_tclust_fixed(x,k,alpha,c,TSR_V2$nstart,TSR_CONFIG))
  if (!macro$ok) stop(name, " TCLUST failed: ", macro$error)
  raw <- as.integer(macro$value$cluster)

  set.seed(TSR_V2$base_seed + 200L)
  micro <- run_safe(fit_tsr_from_macro(x,raw,k,validated=TRUE,null_reps=TSR_V2$null_reps,
                                       structure_alpha=TSR_V2$structure_alpha,
                                       min_density_ratio=TSR_V2$min_density_ratio,
                                       novelty_prob=TSR_V2$novelty_prob,
                                       eps_quantile=TSR_V2$eps_quantile,
                                       eps_grid_size=TSR_V2$eps_grid_size))
  if (!micro$ok) stop(name, " TSRClust failed: ", micro$error)
  final <- micro$value$cluster

  set.seed(TSR_V2$base_seed + 201L)
  micro_unval <- run_safe(
    fit_tsr_from_macro(
      x, raw, k,
      validated = FALSE,
      null_reps = TSR_V2$null_reps,
      structure_alpha = TSR_V2$structure_alpha,
      min_density_ratio = TSR_V2$min_density_ratio,
      novelty_prob = TSR_V2$novelty_prob,
      eps_quantile = TSR_V2$eps_quantile,
      eps_grid_size = TSR_V2$eps_grid_size
    )
  )

  if (!micro_unval$ok)
    stop(name, " TSRClust unvalidated failed: ", micro_unval$error)

  final_unval <- micro_unval$value$cluster
  r0 <- which(raw==0L)

  timings <- data.frame(
    dataset=name,
    method=c("TCLUST","TSRClust_validated"),
    preprocess_seconds=prep,
    macro_seconds=c(macro$elapsed,macro$elapsed),
    micro_seconds=c(0,micro$elapsed),
    end_to_end_seconds=c(prep+macro$elapsed,prep+macro$elapsed+micro$elapsed),
    stringsAsFactors=FALSE
  )

  # Full-data DBSCAN/HDBSCAN use the same standardized input; package loading and output writing are excluded.
  timings <- rbind(
    timings,
    data.frame(
      dataset = name,
      method = "TSRClust_unvalidated",
      preprocess_seconds = prep,
      macro_seconds = macro$elapsed,
      micro_seconds = micro_unval$elapsed,
      end_to_end_seconds = prep + macro$elapsed + micro_unval$elapsed,
      stringsAsFactors = FALSE
    )
  )

  db <- run_safe(fit_full_dbscan(x, TSR_V2$eps_quantile))
  if (identical(name, "Statlog_Shuttle")) {
    hd <- list(
      ok = FALSE,
      elapsed = NA_real_,
      error = "Skipped: HDBSCAN caused segmentation fault on full Shuttle data.",
      warning = character()
    )
  } else {
    hd <- run_safe(fit_full_hdbscan(x))
  }
  if (db$ok) timings <- rbind(timings,data.frame(dataset=name,method="DBSCAN_full",preprocess_seconds=prep,
       macro_seconds=0,micro_seconds=0,end_to_end_seconds=prep+db$elapsed))
  if (hd$ok) timings <- rbind(timings,data.frame(dataset=name,method="HDBSCAN_full",preprocess_seconds=prep,
       macro_seconds=0,micro_seconds=0,end_to_end_seconds=prep+hd$elapsed))

  method_labels <- list(
    TCLUST = raw,
    TSRClust_validated = final,
    TSRClust_unvalidated = final_unval
  )
  if (db$ok) method_labels$DBSCAN_full <- db$value$cluster
  if (hd$ok) method_labels$HDBSCAN_full <- hd$value$cluster

  if (run_global_robust) {
    set.seed(TSR_V2$base_seed + 300L)
    oo <- run_safe(fit_otrimle_fixed(x,G=k,erc=TSR_V2$otrimle_erc,npr.max=TSR_V2$otrimle_npr_max))
    if (oo$ok) {
      method_labels$OTRIMLE_fixedK <- oo$value$cluster
      timings <- rbind(timings,data.frame(dataset=name,method="OTRIMLE_fixedK",preprocess_seconds=prep,
        macro_seconds=0,micro_seconds=0,end_to_end_seconds=prep+oo$elapsed))
    }
    set.seed(TSR_V2$base_seed + 400L)
    tg <- run_safe(fit_teigen_bic(x,Gs=1:min(6L,k+3L)))
    if (tg$ok) {
      method_labels$tEIGEN_BIC <- tg$value$cluster
      timings <- rbind(timings,data.frame(dataset=name,method="tEIGEN_BIC",preprocess_seconds=prep,
        macro_seconds=0,micro_seconds=0,end_to_end_seconds=prep+tg$elapsed))
    }
  }

  macro_ari <- do.call(rbind,lapply(names(method_labels),function(m){
    data.frame(dataset=name,method=m,MacroLabelARI=simple_ari(macro_truth,method_labels[[m]]),stringsAsFactors=FALSE)
  }))

  residual_eval <- do.call(rbind,lapply(names(method_labels),function(m){
    lab <- method_labels[[m]][r0]
    data.frame(dataset=name,method=m,n_residual=length(r0),
      ResidualLocalARI=if(length(r0)>1) simple_ari(local_truth[r0],lab) else NA_real_,
      ResidualAssignedFraction=if(length(r0)) mean(lab!=0L) else NA_real_,stringsAsFactors=FALSE)
  }))

  trim_region <- data.frame(dataset=name,local_class=levels(factor(local_truth)),stringsAsFactors=FALSE)
  trim_region$total <- as.integer(table(factor(local_truth,levels=trim_region$local_class)))
  trim_region$trimmed <- as.integer(table(factor(local_truth[r0],levels=trim_region$local_class)))
  trim_region$trim_rate <- trim_region$trimmed/trim_region$total

  micro_comp <- if(length(r0)) cluster_purity(local_truth[r0], final[r0], positive_only=TRUE) else data.frame()
  if(nrow(micro_comp)) micro_comp$dataset <- name

  list(x=x,raw=raw,final=final,r0=r0,timings=timings,macro_ari=macro_ari,
       residual_eval=residual_eval,trim_region=trim_region,micro_comp=micro_comp,
       micro_info=micro$value$micro)
}

# 1) Italian olive-oil chemometric data.
data("oliveoil", package="pdfCluster", envir=environment())
olive_macro <- oliveoil[[1]]
olive_region <- oliveoil[[2]]
olive_raw <- as.matrix(oliveoil[,3:10])
olive_ilr <- ilr_transform(olive_raw, pseudocount=1)
olive_x_std <- robust_standardize(olive_ilr)
olive_c <- TSR_V2$olive_c
olive <- analyse_one("Italian_Olive_Oil", olive_ilr, olive_macro, olive_region,
                     TSR_V2$olive_k, TSR_V2$olive_alpha, olive_c, run_global_robust=TRUE)

# 2) Statlog Shuttle. K=2, alpha=.05, c=20 are frozen from the label-free V2.1 audit.
data("Shuttle", package="mlbench", envir=environment())
shuttle_x <- as.matrix(Shuttle[,1:9])
shuttle_class <- Shuttle[[10]]
# Local and macro external labels are both the 7 observed operational labels here; they are never used in fitting.
shuttle <- analyse_one("Statlog_Shuttle", shuttle_x, shuttle_class, shuttle_class,
                       TSR_V2$shuttle_k, TSR_V2$shuttle_alpha, TSR_V2$shuttle_c,
                       run_global_robust=FALSE)

write.csv(rbind(olive$timings,shuttle$timings),"v2/results/realdata/timings.csv",row.names=FALSE)
write.csv(rbind(olive$macro_ari,shuttle$macro_ari),"v2/results/realdata/external_macro_ari.csv",row.names=FALSE)
write.csv(rbind(olive$residual_eval,shuttle$residual_eval),"v2/results/realdata/residual_external_eval.csv",row.names=FALSE)
write.csv(rbind(olive$trim_region,shuttle$trim_region),"v2/results/realdata/class_trimming_rates.csv",row.names=FALSE)
if(nrow(olive$micro_comp)) write.csv(olive$micro_comp,"v2/results/realdata/olive_microcluster_composition.csv",row.names=FALSE)
if(nrow(shuttle$micro_comp)) write.csv(shuttle$micro_comp,"v2/results/realdata/shuttle_microcluster_composition.csv",row.names=FALSE)

info <- data.frame(dataset=c("Italian_Olive_Oil","Statlog_Shuttle"),
                   n=c(nrow(olive_ilr),nrow(shuttle_x)),p=c(ncol(olive_ilr),ncol(shuttle_x)),
                   k=c(TSR_V2$olive_k,TSR_V2$shuttle_k),alpha=c(TSR_V2$olive_alpha,TSR_V2$shuttle_alpha),
                   c=c(olive_c,TSR_V2$shuttle_c),residual_n=c(length(olive$r0),length(shuttle$r0)),
                   micro_q=c(length(unique(olive$final[olive$final>TSR_V2$olive_k])),
                             length(unique(shuttle$final[shuttle$final>TSR_V2$shuttle_k]))),
                   micro_p_value=c(olive$micro_info$p_value,shuttle$micro_info$p_value))
write.csv(info,"v2/results/realdata/dataset_fit_summary.csv",row.names=FALSE)
print(info,row.names=FALSE)
cat("REALDATA_V2_OK\n")
