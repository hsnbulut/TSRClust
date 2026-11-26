# Reviewer 1 Comment 1 only. Read-only source audit; no fitting or simulation.
stopifnot(as.character(getRversion())=='4.5.2')
filearg<-grep('^--file=',commandArgs(),value=TRUE)
code<-dirname(normalizePath(sub('^--file=','',filearg[1])))
comment<-dirname(code);revision<-dirname(dirname(comment));project<-dirname(revision)
evidence<-file.path(comment,'evidence');tables<-file.path(comment,'tables')
input<-read.csv(file.path(evidence,'source_inventory.csv'),stringsAsFactors=FALSE)
sha256<-function(path) {
 exe<-Sys.which('shasum');args<-c('-a','256',shQuote(path))
 if(!nzchar(exe)){exe<-Sys.which('sha256sum');args<-shQuote(path)}
 if(!nzchar(exe))stop('A read-only SHA256 utility is required')
 strsplit(system2(exe,args,stdout=TRUE),'[[:space:]]+')[[1]][1]
}
resolved<-ifelse(nzchar(input$project_relative),file.path(project,input$project_relative),input$path)
archive_override<-Sys.getenv('TSR_TIER1_ARCHIVE')
if(nzchar(archive_override))resolved[input$role=='immutable Tier1 archive']<-archive_override
for(i in seq_len(nrow(input)))stopifnot(file.exists(resolved[i]),sha256(resolved[i])==input$sha256[i])
archive<-resolved[input$role=='immutable Tier1 archive'];stopifnot(length(archive)==1L)
work<-tempfile('comment1_tier1_');dir.create(work)
main<-function() {
 on.exit(unlink(work,recursive=TRUE),add=TRUE)
 utils::untar(archive,files=c('replications.csv','summary.csv','paired_differences.csv','micro_counts.csv'),exdir=work)
 x<-read.csv(file.path(work,'replications.csv'),stringsAsFactors=FALSE,check.names=FALSE)
 summary<-read.csv(file.path(work,'summary.csv'),stringsAsFactors=FALSE)
 design<-read.csv(file.path(project,'revision_scientific_reports/config/design.csv'),stringsAsFactors=FALSE);design<-design[design$tier==1,]
 tasks<-read.csv(file.path(project,'revision_scientific_reports/config/tasks_tier1.csv'))
 stopifnot(nrow(x)==98800L,nrow(design)==282L,nrow(tasks)==9880L,all(x$status=='OK'),!anyDuplicated(paste(x$cell_id,x$replication)))
 expected<-unlist(lapply(seq_len(nrow(tasks)),function(i)paste(tasks$cell_id[i],seq.int(tasks$rep_start[i],tasks$rep_end[i]))))
 stopifnot(setequal(expected,paste(x$cell_id,x$replication)),setequal(x$cell_id,design$cell_id))
 for(nm in names(design)) {
  ii<-match(x$cell_id,design$cell_id);stopifnot(identical(as.character(x[[nm]]),as.character(design[[nm]][ii])))
 }
 seed<-function(dataset,rep,stream,repeat_id=1)stream*600000000+(repeat_id-1)*60000000+dataset*10000+rep
 stopifnot(all(x$data_seed==seed(x$dataset_id,x$replication,0)),all(x$algorithm_seed==seed(x$dataset_id,x$replication,1,x$algorithm_repeat)),all(x$null_seed==seed(x$dataset_id,x$replication,2,x$null_repeat)))
 leg<-new.env(parent=globalenv())
 sys.source(file.path(project,'Scientific Reports Manuscript/reproducibility/core/utils.R'),envir=leg)
 sys.source(file.path(project,'Scientific Reports Manuscript/reproducibility/evaluation/evaluation_metrics.R'),envir=leg)
 ext<-new.env(parent=globalenv());ext$legacy<-leg
 ext$check_labels<-function(truth,pred)stopifnot(length(truth)==length(pred),!anyNA(pred),all(pred>=0))
 sys.source(file.path(project,'revision_scientific_reports/R/metrics.R'),envir=ext)
 metrics<-names(leg$evaluate_partition(c(1,1,2,2,3,3,0),c(1,1,2,2,3,3,0),2))
 additional<-c('MinorityClusterF1_legacy_v1','MinorityDetectionF1_legacy_v1','MinorityPrecision_legacy_v1','MinorityPrecision_all_extra_v2','MinorityRecall_all_extra_v2','MinorityDetectionF1_all_extra_v2','MacroFalsePositiveReassignment','any_native_micro','n_accepted_micro','residual_assignment_rate')
 metrics<-c(metrics,additional);stopifnot(all(metrics%in%names(x)))
 available<-data.frame(metric=metrics,n_rows=nrow(x),n_finite=vapply(x[metrics],function(z)sum(is.finite(z)),integer(1)),n_undefined=vapply(x[metrics],function(z)sum(!is.finite(z)),integer(1)))
 available$n_macro_applicable<-ifelse(available$metric%in%c('MajorARI','MajorClusterF1','MajorRetention','MacroFalsePositiveReassignment'),sum(x$mode!='residual_only'),NA_integer_)
 write.csv(available,file.path(tables,'tier1_metric_availability.csv'),row.names=FALSE)
 families<-aggregate(list(n_rows=x$cell_id),x[c('experiment','method','mode','regime')],length)
 write.csv(families,file.path(tables,'tier1_method_experiment_counts.csv'),row.names=FALSE)
 write.csv(as.data.frame(table(status=x$status)),file.path(evidence,'tier1_status_counts.csv'),row.names=FALSE)
 write.csv(data.frame(column=names(x),storage_type=vapply(x,class,character(1))),file.path(evidence,'tier1_columns.csv'),row.names=FALSE)
 checks<-data.frame(check=character(),status=character(),n=integer(),detail=character())
 pass<-function(name,n,detail=''){checks<<-rbind(checks,data.frame(check=name,status='PASS',n=n,detail=detail));cat('PASS',name,n,'\n')}
 pass('complete unique cell replication coverage',nrow(x),'9880 original ten-fit blocks; 282 cells; all expected rows present once')
 pass('frozen design and deterministic seeds',nrow(x),'All design columns and three seeds agree with frozen source')
 eq<-function(a,b)all((is.na(a)&is.na(b))|(!is.na(a)&!is.na(b)&abs(a-b)<1e-10))
 for(pair in list(c('MinorityClusterF1','MinorityClusterF1_legacy_v1'),c('MinorityDetectionF1','MinorityDetectionF1_legacy_v1'),c('MinorityPrecision','MinorityPrecision_legacy_v1'))) {
  stopifnot(eq(x[[pair[1]]],x[[pair[2]]]));pass(paste(pair,collapse=' equals '),nrow(x))
 }
 stopifnot(all(x$AnyMicroDiscovery==as.numeric(x$QHat>0)),all(is.na(x$MinorityDetectionF1[x$QTrue==0])),all(is.na(x$MinorityDetectionF1_all_extra_v2[x$QTrue==0])),all(is.na(x$MinorityClusterF1[x$QTrue==0])))
 pass('absence of a true minority leaves minority F1 undefined',sum(x$QTrue==0),'Undefined values are retained, not replaced by zero')
 positive<-x$QTrue>0;finite<-positive&is.finite(x$MinorityDetectionF1)&is.finite(x$MinorityDetectionF1_all_extra_v2)
 stopifnot(all(x$MinorityRecall_all_extra_v2[positive]+1e-10>=x$MinorityRecall[positive]))
 audit<-do.call(rbind,lapply(split(which(finite),x$method[finite]),function(ix){
  delta<-x$MinorityDetectionF1_all_extra_v2[ix]-x$MinorityDetectionF1[ix]
  data.frame(method=x$method[ix[1]],n_positive=length(ix),n_different=sum(abs(delta)>1e-10),n_v2_lower=sum(delta< -1e-10),n_v2_higher=sum(delta>1e-10),delta_min=min(delta),delta_max=max(delta),legacy_mean=mean(x$MinorityDetectionF1[ix]),all_extra_mean=mean(x$MinorityDetectionF1_all_extra_v2[ix]))
 }))
 write.csv(audit,file.path(tables,'detection_metric_version_audit.csv'),row.names=FALSE)
 native<-x$method=='TSRClust'&is.finite(x$any_native_micro)
 agreement<-aggregate(list(n=native[native],different=x$AnyMicroDiscovery[native]!=x$any_native_micro[native]),x[native,c('mode','regime')],sum)
 write.csv(agreement,file.path(tables,'discovery_definition_audit.csv'),row.names=FALSE)
 pass('AnyMicroDiscovery equals eligible extra count event',nrow(x),'Native Stage-II discovery separately audited; not assumed equivalent')
 # Validate every exported scalar summary against authoritative replication rows.
 err<-numeric(nrow(summary));nfinite<-integer(nrow(summary))
 for(i in seq_len(nrow(summary))) {
  v<-x[x$cell_id==summary$cell_id[i],summary$metric[i]];v<-v[is.finite(v)];nfinite[i]<-length(v)
  got<-if(length(v))mean(v) else NA_real_;expected_mean<-summary$mean[i]
  err[i]<-if(is.na(got)&&is.na(expected_mean))0 else if(is.na(got)||is.na(expected_mean))Inf else abs(got-expected_mean)
 }
 stopifnot(all(err<1e-10),all(nfinite==summary$n_finite),all(summary$n_failed==0),all(summary$n_missing==0))
 pass('all exported means and finite denominators reproduced',nrow(summary),paste('maximum absolute difference',format(max(err),digits=4)))
 write.csv(data.frame(cell_id=summary$cell_id,metric=summary$metric,absolute_mean_difference=err,n_finite_recomputed=nfinite,n_finite_exported=summary$n_finite),file.path(evidence,'summary_reproduction_checks.csv'),row.names=FALSE)
 # Verify the companion exports, including paired means, without inferring labels.
 micro<-read.csv(file.path(work,'micro_counts.csv'),stringsAsFactors=FALSE)
 jj<-match(paste(x$cell_id,x$replication),paste(micro$cell_id,micro$replication))
 stopifnot(nrow(micro)==nrow(x),!anyNA(jj),all(micro$status[jj]==x$status),eq(micro$n_accepted_micro[jj],x$n_accepted_micro))
 pass('micro-count export matches every replication',nrow(micro))
 paired<-read.csv(file.path(work,'paired_differences.csv'),stringsAsFactors=FALSE)
 perr<-numeric(nrow(paired))
 for(i in seq_len(nrow(paired))) {
  z<-x[x$cell_id==paired$cell_id[i],c('replication',paired$metric[i])]
  b<-x[x$cell_id==paired$baseline_cell[i],c('replication',paired$metric[i])]
  w<-merge(z,b,by='replication');v<-w[[2]]-w[[3]];v<-v[is.finite(v)]
  val<-if(length(v))mean(v) else NA_real_
  perr[i]<-if(is.na(val)&&is.na(paired$mean_difference[i]))0 else abs(val-paired$mean_difference[i])
  stopifnot(length(v)==paired$n_pairs[i])
 }
 stopifnot(all(perr<1e-10));pass('all exported paired means and finite-pair denominators reproduced',nrow(paired),paste('maximum absolute difference',format(max(perr),digits=4)))
 write.csv(data.frame(cell_id=paired$cell_id,baseline_cell=paired$baseline_cell,metric=paired$metric,absolute_mean_difference=perr,n_pairs=paired$n_pairs),file.path(evidence,'paired_reproduction_checks.csv'),row.names=FALSE)
 # Deterministic metric counterexamples, not simulation replications or fitted results.
 truth<-c(rep(1,10),rep(2,10),rep(3,4),rep(0,4))
 examples<-list(unmatched_noise_cluster=list(y=truth,z=c(rep(10,10),rep(20,10),rep(30,4),rep(40,4)),k=2),
  split_minority=list(y=truth,z=c(rep(10,10),rep(20,10),30,30,31,31,rep(0,4)),k=2),
  trimmed_major_class=list(y=c(rep(1,4),rep(2,4)),z=c(rep(0,4),rep(10,4)),k=2))
 toy<-do.call(rbind,lapply(names(examples),function(nm){q<-examples[[nm]];v<-ext$revision_metrics(q$y,q$z,q$k)$metrics;data.frame(example=nm,v,check.names=FALSE)}))
 stopifnot(abs(toy$MinorityDetectionF1[1]-1)<1e-10,abs(toy$MinorityDetectionF1_all_extra_v2[1]-2/3)<1e-10,
  abs(toy$MinorityClusterF1[2]-2/3)<1e-10,abs(toy$MinorityDetectionF1_all_extra_v2[2]-1)<1e-10,
  abs(toy$MajorARI[3]-1)<1e-10,abs(toy$MajorRetention[3]-.5)<1e-10)
 write.csv(toy,file.path(evidence,'definition_counterexamples.csv'),row.names=FALSE)
 pass('source evaluator counterexamples',nrow(toy),'Synthetic fixed label vectors only; no fitting/random generation')
 pkg<-new.env(parent=globalenv())
 sys.source(file.path(project,'Scientific Reports Manuscript/RPackageNew/TSRClust/R/helpers.R'),envir=pkg)
 sys.source(file.path(project,'Scientific Reports Manuscript/RPackageNew/TSRClust/R/metrics.R'),envir=pkg)
 native_labels<-c(rep(1,10),rep(2,10),rep(3,4),rep(4,4))
 package_score<-pkg$evaluate_tsrclust(truth,native_labels,2)$MinorityF1
 workflow<-ext$revision_metrics(truth,native_labels,2)$metrics
 stopifnot(abs(package_score-2/3)<1e-10,abs(workflow$MinorityDetectionF1-1)<1e-10,abs(workflow$MinorityDetectionF1_all_extra_v2-package_score)<1e-10)
 write.csv(data.frame(example='fixed native labels with unmatched extra noise cluster',package_MinorityF1=package_score,workflow_matched_union_F1=workflow$MinorityDetectionF1,workflow_all_extra_v2_F1=workflow$MinorityDetectionF1_all_extra_v2),file.path(evidence,'package_helper_scope_check.csv'),row.names=FALSE)
 pass('package convenience helper is not the cross-method matched-union evaluator',1,'Native-label convention and hierarchical matching are explicitly distinguished')

 # Reuse submitted real-data CSVs; do not refit or reconstruct memberships.
 real<-file.path(project,'Scientific Reports Manuscript/github_results/real_data_results')
 fit<-read.csv(file.path(real,'dataset_fit_summary.csv'));extari<-read.csv(file.path(real,'external_macro_ari.csv'));res<-read.csv(file.path(real,'residual_external_eval.csv'))
 compo<-read.csv(file.path(real,'olive_microcluster_composition.csv'));comps<-read.csv(file.path(real,'shuttle_microcluster_composition.csv'))
 for(i in seq_len(nrow(fit))) {
  comp<-if(fit$dataset[i]=='Italian_Olive_Oil')compo else comps
  idx<-res$dataset==fit$dataset[i]&res$method=='TSRClust_validated'
  stopifnot(nrow(comp)==fit$micro_q[i],abs(sum(comp$size)/fit$residual_n[i]-res$ResidualAssignedFraction[idx])<1e-12)
 }
 real_audit<-merge(extari,res,by=c('dataset','method'));real_audit<-merge(real_audit,fit[c('dataset','n')],by='dataset')
 names(real_audit)[names(real_audit)=='MacroLabelARI']<-'ExternalLabelARI_all_observations'
 real_audit$external_reference<-ifelse(real_audit$dataset=='Italian_Olive_Oil','three macro-areas on all 572 oils','seven operational classes on all 58000 observations')
 real_audit$residual_reference<-ifelse(real_audit$dataset=='Italian_Olive_Oil','nine regions restricted to Stage-I TCLUST R0','seven operational classes restricted to Stage-I TCLUST R0')
 write.csv(real_audit,file.path(tables,'submitted_realdata_metric_scope.csv'),row.names=FALSE)
 pass('submitted real-data assignment denominators',nrow(fit),'30/58 for Olive; 511/2900 for Shuttle; no new real-data fit')
 # Scalar overview suitable for a comment-specific manuscript table.
 f1diff<-sum(abs(x$MinorityDetectionF1_all_extra_v2[finite]-x$MinorityDetectionF1[finite])>1e-10)
 overview<-data.frame(item=c('successful_rows','design_cells','task_blocks','method_labels','experiment_families','macro_ARI_finite','global_ARI_finite','true_minority_rows','no_true_minority_rows','detection_versions_differ','native_vs_external_discovery_differ'),
  value=c(nrow(x),length(unique(x$cell_id)),nrow(tasks),length(unique(x$method)),length(unique(x$experiment)),sum(is.finite(x$MajorARI)),sum(is.finite(x$GlobalARI)),sum(positive),sum(!positive),f1diff,sum(x$AnyMicroDiscovery[native]!=x$any_native_micro[native])))
 write.csv(overview,file.path(tables,'comment1_verified_counts.csv'),row.names=FALSE)
 write.csv(checks,file.path(evidence,'validation_checks.csv'),row.names=FALSE)
 writeLines(capture.output(sessionInfo()),file.path(evidence,'R_session_info.txt'))
 saveRDS(list(counts=overview,availability=available,method_family=families,metric_versions=audit,realdata=real_audit),file.path(evidence,'compact_comment1_analysis.rds'))
 cat('COMMENT1_ANALYSIS_OK; source rows',nrow(x),'cells',nrow(design),'minority rows',sum(positive),'version differences',f1diff,'\n')
}
main()
