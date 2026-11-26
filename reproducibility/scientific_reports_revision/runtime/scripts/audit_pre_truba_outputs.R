source('R/load.R')
out<-'results/pilot_R452_pre_truba';z<-read.csv(file.path(out,'replications.csv'),stringsAsFactors=FALSE)
old<-read.csv('results/pilot/replications.csv',stringsAsFactors=FALSE)
failed<-old[old$status=='ERROR',c('cell_id','replication','method','error')]
w<-merge(failed,z[c('cell_id','replication','status','error','data_seed','algorithm_seed','null_seed')],by=c('cell_id','replication'),suffixes=c('_previous','_rerun'))
write.csv(w,'results/pre_truba/previous_eight_reruns.csv',row.names=FALSE)
stopifnot(nrow(w)==8,all(w$status=='OK'))
# Independently check denominators and matching on every fitted partition.
checks<-list();contingencies<-list()
for(f in list.files(out,pattern='^task_.*[.]rds$',full.names=TRUE)) {
 obj<-readRDS(f)
 for(v in obj$results) {
  if(v$row$status!='OK')next
  y<-v$diagnostics$truth;pred<-v$diagnostics$predicted;k<-v$row$main_k
  mt<-legacy$hierarchical_match(y,pred,k)
  for(i in seq_along(mt$true_minority)) {
   g<-mt$matched_minority_pred[i];u<-mt$true_minority[i]
   value<-if(!g)0 else 2*sum(y==u&pred==g)/(sum(y==u)+sum(pred==g))
   saved<-if(!mt$minority_cols[i])0 else mt$minority_score[i,mt$minority_cols[i]]
   stopifnot(abs(value-saved)<1e-12)
  }
  positive<-pred%in%mt$remaining;tp<-sum(positive&y>k);fp<-sum(positive&y<=k);fn<-sum(!positive&y>k)
  f1<-if(any(y>k))legacy$safe_ratio(2*tp,2*tp+fp+fn,0) else NA_real_
  stopifnot(isTRUE(all.equal(f1,v$row$MinorityDetectionF1_all_extra_v2)))
  checks[[length(checks)+1L]]<-data.frame(cell_id=v$row$cell_id,replication=v$row$replication,method=v$row$method,MajorARI=v$row$MajorARI,
   MinorityClusterF1_legacy_v1=v$row$MinorityClusterF1_legacy_v1,MinorityDetectionF1_legacy_v1=v$row$MinorityDetectionF1_legacy_v1,
   MinorityDetectionF1_all_extra_v2=f1,extra_false_positive_observations=fp,macro_observations_in_extras=sum(positive&y%in%seq_len(k)),denominators_verified=TRUE)
  t<-as.data.frame(table(truth=y,predicted=pred));t$cell_id<-v$row$cell_id;t$replication<-v$row$replication;t$method<-v$row$method;contingencies[[length(contingencies)+1L]]<-t
 }
}
write.csv(bind_rows(checks),'results/pre_truba/all_method_metric_integrity.csv',row.names=FALSE)
write.csv(bind_rows(contingencies),'results/pre_truba/all_method_contingencies.csv',row.names=FALSE)
# Repeat one optional dependency fit with the same dataset/algorithm/null seeds.
design<-read.csv('config/design.csv');mcd<-design[design$cell_id==80,];a<-fit_revision(mcd,1);b<-fit_revision(mcd,1)
stopifnot(a$row$status=='OK',b$row$status=='OK',identical(a$diagnostics$predicted,b$diagnostics$predicted),identical(a$diagnostics$micro$null_scores,b$diagnostics$micro$null_scores))
# All task completions are accounted for.
stopifnot(nrow(z)==40,all(z$status=='OK'))
write.csv(data.frame(check=c('all eight formerly failed fits now succeed','all 40 fitted metric denominators checked','MCD repeat deterministic'),status='PASS'),'results/pre_truba/output_checks.csv',row.names=FALSE)
cat('PRE_TRUBA_OUTPUT_CHECKS_OK\n')
