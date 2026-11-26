source('R/load.R')
a<-commandArgs(TRUE);if(length(a)<2)stop('Usage aggregate.R manifest output_dir')
tasks<-read.csv(a[1]);out<-a[2];artifact<-Sys.getenv('TSR_ARTIFACT_DIR');table_dir<-if(nzchar(artifact))file.path(artifact,'tables') else 'tables';figure_dir<-if(nzchar(artifact))file.path(artifact,'figures') else 'figures';dir.create(table_dir,recursive=TRUE,showWarnings=FALSE);dir.create(figure_dir,recursive=TRUE,showWarnings=FALSE);design<-read.csv('config/design.csv');xs<-list();diags<-list()
for(i in seq_len(nrow(tasks))) {
 task<-tasks[i,];path<-file.path(out,sprintf('task_%06d.rds',task$task_id));obj<-if(file.exists(path))readRDS(path) else NULL
 for(rep in seq.int(task$rep_start,task$rep_end)) {
  z<-if(!is.null(obj))obj$results[[as.character(rep)]] else NULL
  if(is.null(z))z<-list(row=cbind(design[design$cell_id==task$cell_id,],data.frame(replication=rep,status='MISSING',error='No result file/replication')))
  xs[[length(xs)+1L]]<-z$row
  diags[[paste(task$cell_id,rep,sep='_')]]<-z$diagnostics
 }
}
raw<-bind_rows(xs);atomic_rds(raw,file.path(out,'replications.rds'));write.csv(raw,file.path(out,'replications.csv'),row.names=FALSE)
metrics<-intersect(c('MajorARI','MinorityClusterF1','MinorityDetectionF1','MinorityPrecision','MinorityRecall',
 'MinorityDetectionF1_all_extra_v2','MinorityPrecision_all_extra_v2','MacroFalsePositiveReassignment','GlobalARI','AbsQError','SpuriousExtraClusters',
 'AnyMicroDiscovery','any_native_micro','n_accepted_micro','residual_assignment_rate','end_to_end_seconds'),names(raw))
summary<-list()
for(id in unique(raw$cell_id)) {
 z<-raw[raw$cell_id==id,];ok<-z$status=='OK'
 for(m in metrics) {
  v<-z[[m]][ok];v<-v[is.finite(v)];nn<-length(v);ci<-if(nn>1)mean(v)+c(-1,1)*qt(.975,nn-1)*sd(v)/sqrt(nn) else c(NA,NA)
  binary<-m %in% c('AnyMicroDiscovery','any_native_micro');events<-if(binary&&nn)sum(v) else NA_real_
  if(binary&&nn)ci<-exact_interval(events,nn)
  summary[[length(summary)+1L]]<-data.frame(cell_id=id,metric=m,n_expected=nrow(z),n_success=sum(ok),n_failed=sum(z$status=='ERROR'),n_missing=sum(z$status=='MISSING'),n_finite=nn,events=events,
   mean=if(nn)mean(v) else NA_real_,sd=if(nn>1)sd(v) else NA_real_,mcse=if(nn>1)sd(v)/sqrt(nn) else NA_real_,lower=ci[1],upper=ci[2],
   median=if(nn)median(v) else NA_real_,q25=if(nn)unname(quantile(v,.25)) else NA_real_,q75=if(nn)unname(quantile(v,.75)) else NA_real_,
   failure_rate_lower=if(binary&&nn)events/nrow(z) else NA_real_,failure_rate_upper=if(binary&&nn)(events+sum(!ok))/nrow(z) else NA_real_)
 }
}
sumtab<-bind_rows(summary);write.csv(sumtab,file.path(out,'summary.csv'),row.names=FALSE)
# Paired parameter effects: match on DGP and independent RNG-repeat IDs, never on scores.
paired<-list()
base<-raw[raw$method=='TSRClust'&raw$rule=='legacy_rule'&raw$reference=='legacy_empirical'&raw$novelty==.995&raw$density==2&raw$minpts_multiplier==2&raw$B==49&raw$regime!='data_selected',]
for(id in unique(raw$cell_id)) {
 z<-raw[raw$cell_id==id,];b<-base[base$dataset_id==z$dataset_id[1]&base$algorithm_repeat==z$algorithm_repeat[1]&base$null_repeat==z$null_repeat[1],]
 if(!nrow(b))next
 b<-b[!duplicated(b$replication),];w<-merge(z,b,by='replication',suffixes=c('_new','_base'))
 for(m in metrics)if(all(paste0(m,c('_new','_base'))%in%names(w))) {
  good<-w$status_new=='OK'&w$status_base=='OK';v<-w[[paste0(m,'_new')]][good]-w[[paste0(m,'_base')]][good];v<-v[is.finite(v)];n<-length(v)
  ci<-if(n>1)mean(v)+c(-1,1)*qt(.975,n-1)*sd(v)/sqrt(n) else c(NA,NA)
  paired[[length(paired)+1L]]<-data.frame(cell_id=id,baseline_cell=if(nrow(b))b$cell_id[1] else NA,metric=m,n_pairs=n,mean_difference=if(n)mean(v) else NA,lower=ci[1],upper=ci[2])
 }
}
write.csv(bind_rows(paired),file.path(out,'paired_differences.csv'),row.names=FALSE)
# Each row is a traceable performance result; TeX is generated without conclusions.
headline<-sumtab[sumtab$metric%in%c('any_native_micro','MinorityClusterF1','MinorityDetectionF1_all_extra_v2','end_to_end_seconds'),]
write.csv(headline,file.path(table_dir,paste0(basename(out),'_summary.csv')),row.names=FALSE)
tex<-c('\\begin{tabular}{rlrrrr}','Cell & Metric & Mean & Lower & Upper & $n$ \\\\','\\hline',
 apply(headline,1,function(z)sprintf('%s & %s & %s & %s & %s & %s \\\\',z['cell_id'],gsub('_','\\_',z['metric'],fixed=TRUE),z['mean'],z['lower'],z['upper'],z['n_finite'])),'\\end{tabular}')
writeLines(tex,file.path(table_dir,paste0(basename(out),'_summary.tex')))
for(metric in c('any_native_micro','MinorityClusterF1','end_to_end_seconds')) {
 z<-sumtab[sumtab$metric==metric&is.finite(sumtab$mean),];if(!nrow(z))next
 pdf(file.path(figure_dir,paste0(basename(out),'_',metric,'.pdf')),width=9,height=5)
 plot(z$cell_id,z$mean,xlab='Design cell ID',ylab=metric,pch=19);good<-is.finite(z$lower)&is.finite(z$upper)&z$upper>z$lower
 arrows(z$cell_id[good],z$lower[good],z$cell_id[good],z$upper[good],angle=90,code=3,length=.03);dev.off()
}
# Explicitly unmeasured cells have NA runtime; extrapolation is only within a measured cell.
work<-read.csv('tables/workload_by_cell.csv');times<-raw[raw$status=='OK'&is.finite(raw$end_to_end_seconds),]
if(nrow(times)) {
 measured<-aggregate(cbind(end_to_end_seconds,result_object_bytes)~cell_id,times,mean)
 w<-merge(work,measured,by='cell_id',all.x=TRUE);w$estimated_seconds<-w$outer_fits*w$end_to_end_seconds;w$estimated_uncompressed_bytes<-w$outer_fits*w$result_object_bytes
 write.csv(w,file.path(table_dir,'workload_measured.csv'),row.names=FALSE)
}
# RNG repeats: paired 2x2 contrasts within datasets, not extra independent datasets.
rng<-raw[raw$experiment=='C-RNG'&raw$status=='OK',]
if(nrow(rng)) {
 rngout<-list()
 for(sc in unique(rng$scenario))for(m in intersect(metrics,names(rng))) {
  z<-rng[rng$scenario==sc,];bydata<-aggregate(z[[m]],list(replication=z$replication),mean,na.rm=TRUE)
  within<-aggregate(z[[m]],list(replication=z$replication),var,na.rm=TRUE)
  if(any(is.finite(bydata$x)))rngout[[length(rngout)+1L]]<-data.frame(scenario=sc,metric=m,n_datasets=sum(is.finite(bydata$x)),between_dataset_variance=var(bydata$x,na.rm=TRUE),mean_within_dataset_rng_variance=mean(within$x,na.rm=TRUE))
 }
 write.csv(bind_rows(rngout),file.path(out,'rng_variance_components_descriptive.csv'),row.names=FALSE)
}
# Raw accepted-count distribution and rejection diagnostics remain usable without large fit objects.
if('n_accepted_micro'%in%names(raw))write.csv(raw[intersect(c('cell_id','replication','status','n_accepted_micro','rejection_reason'),names(raw))],file.path(out,'micro_counts.csv'),row.names=FALSE)
cat('AGGREGATION_OK',nrow(raw),'expected rows,',sum(raw$status=='OK'),'successes\n')
