# Read-only metric/seed verification; writes diagnostics inside this new run only.
source('R/load.R')
a<-commandArgs(TRUE);stopifnot(length(a)==2L)
tasks<-read.csv(a[1]);out<-a[2];expected_identity<-list(manifest_md5=unname(tools::md5sum(a[1])),design_md5=unname(tools::md5sum('config/design.csv')),runtime_md5=unname(tools::md5sum('config/runtime_md5.csv')),version=REVISION_VERSION);design<-read.csv('config/design.csv')
stopifnot(nrow(tasks)==20L,sum(tasks$rep_end-tasks$rep_start+1L)==40L)
checks<-list();cont<-list();missing<-failed<-0L
for(i in seq_len(nrow(tasks))) {
 task<-tasks[i,];path<-file.path(out,sprintf('task_%06d.rds',task$task_id));obj<-if(file.exists(path))readRDS(path) else NULL
 if(!is.null(obj))stopifnot(identical(obj$identity,expected_identity))
 cell<-design[design$cell_id==task$cell_id,,drop=FALSE]
 for(rep in seq.int(task$rep_start,task$rep_end)) {
  z<-if(is.null(obj))NULL else obj$results[[as.character(rep)]]
  if(is.null(z)){missing<-missing+1L;next}
  if(z$row$status!='OK'){failed<-failed+1L;next}
  stopifnot(z$row$data_seed==seed_id(cell$dataset_id,rep,0),z$row$algorithm_seed==seed_id(cell$dataset_id,rep,1,cell$algorithm_repeat),z$row$null_seed==seed_id(cell$dataset_id,rep,2,cell$null_repeat))
  y<-z$diagnostics$truth;pred<-z$diagnostics$predicted;k<-cell$main_k
  mt<-legacy$hierarchical_match(y,pred,k)
  matched<-vapply(seq_along(mt$true_minority),function(j){v<-mt$matched_minority_pred[j];u<-mt$true_minority[j];if(!v)0 else 2*sum(y==u&pred==v)/(sum(y==u)+sum(pred==v))},numeric(1))
  expected_cluster<-if(length(matched))mean(matched) else NA_real_
  positive<-pred%in%mt$remaining;tp<-sum(positive&y>k);fp<-sum(positive&y<=k);fn<-sum(!positive&y>k)
  expected_detection<-if(any(y>k))legacy$safe_ratio(2*tp,2*tp+fp+fn,0) else NA_real_
  stopifnot(isTRUE(all.equal(expected_cluster,z$row$MinorityClusterF1)),isTRUE(all.equal(expected_detection,z$row$MinorityDetectionF1_all_extra_v2)))
  stopifnot(is.finite(z$row$replication_cpu_seconds),is.finite(z$row$replication_elapsed_seconds))
  checks[[length(checks)+1L]]<-data.frame(task_id=task$task_id,cell_id=cell$cell_id,replication=rep,method=cell$method,rule=cell$rule,seeds_verified=TRUE,metric_denominators_verified=TRUE)
  t<-as.data.frame(table(truth=y,predicted=pred));t$task_id<-task$task_id;t$replication<-rep;cont[[length(cont)+1L]]<-t
 }
}
write.csv(bind_rows(checks),file.path(out,'acceptance_metric_seed_checks.csv'),row.names=FALSE)
write.csv(bind_rows(cont),file.path(out,'acceptance_contingencies.csv'),row.names=FALSE)
write.csv(data.frame(expected=40,successful=length(checks),failed=failed,missing=missing),file.path(out,'acceptance_fit_counts.csv'),row.names=FALSE)
if(length(checks)!=40||failed||missing)stop('PILOT_ACCEPTANCE_FAILED: failed/missing fits retained in counts/results')
cat('PILOT_40_FITS_METRICS_SEEDS_OK\n')
