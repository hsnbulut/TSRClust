# Read-only Reviewer 1 Comment 3 audit. No model fitting or package installation.
stopifnot(as.character(getRversion())=='4.5.2')
f<-grep('^--file=',commandArgs(),value=TRUE)[1];code<-dirname(normalizePath(sub('^--file=','',f)))
c3<-dirname(code);rev<-dirname(dirname(c3));project<-dirname(rev)
write<-function(x,n,where='tables')write.csv(x,file.path(c3,where,n),row.names=FALSE,na='NA')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
inv<-read.csv(file.path(c3,'evidence/source_inventory.csv'))
input_path<-function(i)if(nzchar(inv$project_relative[i]))file.path(project,inv$project_relative[i]) else if(basename(inv$path[i])=='TSRClust_Tier1_Analysis.tar.gz')archive else inv$path[i]
stopifnot(all(vapply(seq_len(nrow(inv)),function(i)sha(input_path(i))==inv$sha256[i],logical(1))))
legacy<-new.env();sys.source(file.path(project,'Scientific Reports Manuscript/reproducibility/core/utils.R'),legacy);sys.source(file.path(project,'Scientific Reports Manuscript/reproducibility/evaluation/evaluation_metrics.R'),legacy)
ext<-new.env();ext$legacy<-legacy;ext$check_labels<-function(a,b)stopifnot(length(a)==length(b),!anyNA(b),all(b>=0));sys.source(file.path(project,'revision_scientific_reports/R/metrics.R'),ext)
metrics<-c('MajorARI','GlobalARI','MajorRetention','MinorityPrecision','MinorityRecall','MinorityClusterF1','MinorityDetectionF1','MinorityDetectionF1_all_extra_v2','MacroFalsePositiveReassignment','FalseExtraAssignmentRate','SpuriousExtraClusters','NoiseRecall','MacroToMinority_lower','MacroToMinority_upper','MacroToAnyExtra_lower','MacroToAnyExtra_upper','FDP_matched_when_nonempty')
summary_one<-function(v,strata=rep(1,length(v))) {
 ok<-is.finite(v);v<-v[ok];strata<-strata[ok];n<-length(v)
 if(!n)return(c(n_finite=0,n_strata=0,mean=NA,mcse=NA,df=NA,lo=NA,hi=NA))
 mu<-mean(v);parts<-split(v,strata);nc<-lengths(parts)
 vc<-vapply(parts,function(a)if(length(a)>1)length(a)*var(a)/n^2 else NA_real_,numeric(1))
 se<-if(all(is.finite(vc)))sqrt(sum(vc)) else NA_real_
 df<-if(is.finite(se)&&se>0)sum(vc)^2/sum(vc^2/(nc-1)) else NA_real_
 ci<-if(is.finite(se)&&se==0)c(mu,mu) else if(is.finite(df))mu+c(-1,1)*qt(.975,df)*se else c(NA,NA)
 c(n_finite=n,n_strata=length(parts),mean=mu,mcse=se,df=df,lo=ci[1],hi=ci[2])
}
# Set bounds from partition identities and direct rates, never from aggregate ARI/F1.
augment<-function(x,tier=FALSE) {
 nonret<-1-x$MajorRetention
 stopifnot(all(nonret>=0&nonret<=1,na.rm=TRUE))
 noE<-x$QHat==0|x$FalseExtraAssignmentRate==0
 noV<-x$QHat==x$SpuriousExtraClusters|x$QTrue==0
 upperE<-ifelse(noE,0,nonret);upperV<-ifelse(noV,0,upperE)
 x$MacroToAnyExtra_lower<-0;x$MacroToAnyExtra_upper<-upperE
 x$MacroToMinority_lower<-0;x$MacroToMinority_upper<-upperV
 if(tier) {
  x$MacroToAnyExtra_lower<-x$MacroToAnyExtra_upper<-x$MacroFalsePositiveReassignment
  exactV<-noV|x$SpuriousExtraClusters==0|x$MacroFalsePositiveReassignment==0
  x$MacroToMinority_upper<-ifelse(noV,0,x$MacroFalsePositiveReassignment)
  x$MacroToMinority_lower<-ifelse(exactV,x$MacroToMinority_upper,0)
 }
 # Stored precision is undefined/zero with an empty matched union; do not
 # substitute FDP=1 for the zero-detection convention. Positive precision
 # certifies a nonempty minority-matched union in this positive-edge evaluator.
 x$FDP_matched_when_nonempty<-ifelse(is.finite(x$MinorityPrecision)&x$MinorityPrecision>0,1-x$MinorityPrecision,NA_real_)
 x
}
work<-tempfile('c3_tier1_');dir.create(work)
main<-function(){
 on.exit(unlink(work,recursive=TRUE))
 archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
 stopifnot(sha(archive)=='c1852811d7dd94cf8f7b2dbae66d6c17f4284a6cdcf6abbfd8721d73076f94f6')
 untar(archive,files='replications.csv',exdir=work);all<-read.csv(file.path(work,'replications.csv'))
 stopifnot(nrow(all)==98800,all(all$status=='OK'))
 o<-all[all$method=='OTRIMLE_fixedK',];stopifnot(nrow(o)==400)
 cells<-unique(o$cell_id)
 tier<-all[all$method%in%c('OTRIMLE_fixedK','TSRClust','TCLUST')&all$dataset_id%in%unique(o$dataset_id)&all$experiment=='D1-E',]
 stopifnot(nrow(tier)==1200);tier<-augment(tier,TRUE)
 tier$cohort<-paste0('Tier1 ',tier$scenario);tier$source_family<-'Tier1 revision';tier$data_key<-paste(tier$dataset_id,tier$replication)
 original<-list();validation<-list();availability<-list();rawdir<-file.path(project,'Scientific Reports Manuscript/analysis/v2/results')
 for(st in c('architecture','stress','negative_control')) {
  path<-file.path(rawdir,'final_sim_combined',paste0(st,'_replications.rds'));a<-readRDS(path)
  csvpath<-file.path(project,'Scientific Reports Manuscript/github_results/raw_results',paste0(st,'_replications.csv.gz'));b<-read.csv(gzfile(csvpath))
  stopifnot(identical(names(a),names(b)),nrow(a)==nrow(b),all(a$method==b$method),all(a$cell_id==b$cell_id),all(a$replication==b$replication))
  for(m in intersect(metrics,names(a)))stopifnot(identical(is.na(a[[m]]),is.na(b[[m]])),all(abs(a[[m]]-b[[m]])<1e-11,na.rm=TRUE))
  availability[[length(availability)+1]]<-data.frame(path=path,kind='final original RDS and authoritative gzip CSV',rows=nrow(a),columns=ncol(a),list_columns=sum(vapply(a,is.list,logical(1))),observation_labels=FALSE,confusion_counts=FALSE)
  z<-a[a$method%in%c('OTRIMLE_fixedK','TSRClust_validated','TCLUST'),];z$cohort<-if(st=='architecture')ifelse(z$q>0,'Original architecture M1-M4','Original architecture N0') else paste('Original',st)
  z$source_family<-'Original submitted';z$data_key<-paste(st,z$cell_id,z$replication)
  original[[st]]<-augment(z)
 }
 # Inspect every locally available auxiliary worker RDS rather than treating
 # the word 'raw' in a filename as evidence of saved memberships.
 for(folder in c('competitor_pilot_raw','otrimle_sensitivity_raw'))for(path in list.files(file.path(rawdir,folder),pattern='[.]rds$',full.names=TRUE)) {
  a<-readRDS(path);stopifnot(is.data.frame(a),!any(vapply(a,is.list,logical(1))))
  availability[[length(availability)+1]]<-data.frame(path=path,kind=folder,rows=nrow(a),columns=ncol(a),list_columns=0,observation_labels=FALSE,confusion_counts=FALSE)
 }
 write(do.call(rbind,availability),'saved_object_availability.csv','evidence')
 orig<-do.call(rbind,original);orig$MinorityDetectionF1_all_extra_v2<-NA_real_;orig$MacroFalsePositiveReassignment<-NA_real_
 # This NA means unavailable in the ORIGINAL export, not a zero score.
 summaries<-list();scenario<-list();paired<-list();cellpaired<-list();counts<-list();j<-0;k<-0;l<-0;cp<-0
 sources<-list(original=orig,tier1=tier)
 for(sname in names(sources)) {
  x<-sources[[sname]]
  for(cohort in unique(x$cohort))for(method in unique(x$method)) {
   a<-x[x$cohort==cohort&x$method==method,]
   counts[[length(counts)+1]]<-data.frame(source_family=x$source_family[1],cohort=cohort,method=method,n_records=nrow(a),n_OK=sum(a$status=='OK'))
   for(m in metrics){j<-j+1;summaries[[j]]<-cbind(source_family=x$source_family[1],cohort=cohort,method=method,metric=m,n_records=nrow(a),as.data.frame(as.list(summary_one(a[[m]],a$cell_id))))}
  }
  # Per-cell definitions retained; no pooling across original and revision data.
  for(cid in unique(x$cell_id))for(method in unique(x$method)) {
   a<-x[x$cell_id==cid&x$method==method,]
   # Original cell IDs repeat across studies, so split again by cohort.
   if(!nrow(a))next
   for(cohort in unique(a$cohort)) {
    z<-a[a$cohort==cohort,]
    for(m in metrics){k<-k+1;scenario[[k]]<-cbind(source_family=x$source_family[1],cohort=cohort,cell_id=cid,scenario=if('scenario_code'%in%names(z))z$scenario_code[1] else z$scenario[1],factor_changed=if('factor_changed'%in%names(z))z$factor_changed[1] else '',p=z$p[1],separation=z$separation[1],minority_proximity=z$minority_proximity[1],method=method,metric=m,n_records=nrow(z),as.data.frame(as.list(summary_one(z[[m]]))))}
   }
  }
  for(cohort in unique(x$cohort))for(comp in if(sname=='original')c('TSRClust_validated','TCLUST') else c('TSRClust','TCLUST')) {
   a<-x[x$cohort==cohort&x$method=='OTRIMLE_fixedK',];b<-x[x$cohort==cohort&x$method==comp,]
   z<-merge(a,b,by='data_key',suffixes=c('_OT','_comparator'));stopifnot(nrow(z)==nrow(a),!anyDuplicated(z$data_key))
   if(sname=='original')stopifnot(all(z$seed_OT==z$seed_comparator)) else for(ss in c('data_seed','algorithm_seed','null_seed'))stopifnot(all(z[[paste0(ss,'_OT')]]==z[[paste0(ss,'_comparator')]]))
   for(m in setdiff(metrics,c('MacroToMinority_lower','MacroToMinority_upper','MacroToAnyExtra_lower','MacroToAnyExtra_upper'))){
    v<-z[[paste0(m,'_OT')]]-z[[paste0(m,'_comparator')]];l<-l+1
    paired[[l]]<-cbind(source_family=x$source_family[1],cohort=cohort,comparator=comp,direction='OTRIMLE minus comparator',metric=m,n_joined=nrow(z),as.data.frame(as.list(summary_one(v,z$cell_id_OT))),interval='Pointwise stratified Monte Carlo 95%; fixed-cell weighting; Welch degrees of freedom; not multiplicity adjusted; zero variance is not a population guarantee')
    for(cid in unique(z$cell_id_OT)){
     ok<-z$cell_id_OT==cid;cp<-cp+1
     cellpaired[[cp]]<-cbind(source_family=x$source_family[1],cohort=cohort,cell_id=cid,comparator=comp,metric=m,direction='OTRIMLE minus comparator',n_joined=sum(ok),as.data.frame(as.list(summary_one(v[ok]))))
    }
   }
  }
 }
 sumtab<-do.call(rbind,summaries);scen<-do.call(rbind,scenario);pairtab<-do.call(rbind,paired)
 write(sumtab,'method_metric_summaries.csv');write(scen,'scenario_metric_summaries.csv');write(pairtab,'paired_comparisons.csv');write(do.call(rbind,cellpaired),'cell_paired_comparisons.csv');write(do.call(rbind,counts),'record_counts.csv')
 # Reproduce exactly the original reported pooling, without complete-case deletion
 # across unrelated metrics (e.g., NoiseRecall may be undefined without noise).
 reproduction<-list()
 for(st in c('architecture','stress')) {
  a<-original[[st]];a<-a[a$q>0,];pooled<-read.csv(file.path(project,'Scientific Reports Manuscript/github_results/summary_results',if(st=='architecture')'architecture_M1_M4_pooled.csv' else 'stress_pooled.csv'))
  for(method in unique(a$method))for(m in c('MajorARI','MinorityClusterF1','MinorityDetectionF1','SpuriousExtraClusters')) {
   got<-mean(a[[m]][a$method==method],na.rm=TRUE);expected<-pooled[[m]][pooled$method==method]
   stopifnot(length(expected)==1,abs(got-expected)<1e-11)
   reproduction[[length(reproduction)+1]]<-data.frame(study=st,method=method,metric=m,recomputed=got,exported=expected,absolute_difference=abs(got-expected))
  }
 }
 write(do.call(rbind,reproduction),'headline_reproduction_checks.csv','evidence')
 # Diagnosability and interval bounds per source/cohort, never invented point counts.
 diagnostic<-list();j<-0
 for(x in sources)for(cohort in unique(x$cohort)) {
  a<-x[x$cohort==cohort&x$method=='OTRIMLE_fixedK',]
  if(!nrow(a))next
  for(target in c('MacroToMinority','MacroToAnyExtra')) {
   lower<-a[[paste0(target,'_lower')]];upper<-a[[paste0(target,'_upper')]];exact<-is.finite(lower)&is.finite(upper)&lower==upper
   j<-j+1;diagnostic[[j]]<-data.frame(source_family=x$source_family[1],cohort=cohort,target=target,denominator='All true macro observations per fit; summary averages fit-specific proportions',n_records=nrow(a),n_exactly_identified=sum(exact),n_NOT_IDENTIFIABLE=sum(!exact),mean_lower_bound=mean(lower),mean_upper_bound=mean(upper),exact_point_rate_if_all_identified=if(all(exact))mean(upper) else NA_real_,status=if(all(exact))'IDENTIFIED FROM STORED RATES/COUNTS' else 'NOT IDENTIFIABLE: bounds only',basis=if(x$source_family[1]=='Tier1 revision')'Stored MacroFalsePositiveReassignment; V subset E; QHat/Spurious counts. No ARI/F1 inversion.' else '0 <= macro-to-V <= macro-to-E <= 1-MajorRetention; empty sets/zero stored false-extra rate tighten bounds. No ARI/F1 inversion.')
  }
 }
 write(do.call(rbind,diagnostic),'retention_only_promotion_bounds.csv')
 joint<-list()
 for(st in c('architecture','stress')) {
  a<-original[[st]];a<-a[a$method=='OTRIMLE_fixedK'&a$q>0,];z<-a[a$MajorARI<=.1&a$MinorityClusterF1>=.7,]
  joint[[st]]<-data.frame(study=st,rule='Descriptive flag: MajorARI<=0.10 and MinorityClusterF1>=0.70; not a selected validation subset',n_positive=nrow(a),n_flagged=nrow(z),fraction_flagged=nrow(z)/nrow(a),mean_retention_flagged=mean(z$MajorRetention),mean_retention_only_promotion_upper_flagged=mean(z$MacroToMinority_upper))
 }
 write(do.call(rbind,joint),'high_minority_low_macro.csv')
 codes<-do.call(rbind,lapply(names(original),function(st){a<-original[[st]];a<-a[a$method=='OTRIMLE_fixedK',];z<-as.data.frame(table(status=a$status,code=a$otrimle_code,useNA='ifany'));z$study<-st;z}));write(codes,'otrimle_convergence_codes.csv')
 # Fixed label counterexamples, not empirical simulations/fitted data.
 y<-c(rep(1L,10),rep(2L,10),rep(3L,4),rep(0L,4));z<-c(rep(1L,20),rep(2L,4),rep(0L,4))
 contaminated<-c(rep(1L,10),rep(2L,10),rep(3L,4),rep(0L,4));contaminated[1:4]<-3L
 examples<-list(merged_macros_pure_minority=z,contaminated_minority=contaminated,unmatched_noise_extra=replace(z,25:28,3L))
 toy<-do.call(rbind,lapply(names(examples),function(n){p<-examples[[n]];v<-ext$revision_metrics(y,p,2);mt<-v$matching;data.frame(example=n,v$metrics,MacroToMinority=mean(p[y%in%1:2]%in%mt$matched_minority_pred[mt$matched_minority_pred>0]),MacroToAnyExtra=mean(p[y%in%1:2]%in%mt$remaining))}))
 stopifnot(toy$MajorARI[1]==0,toy$MinorityClusterF1[1]==1,toy$MajorRetention[1]==1,toy$MacroToMinority[1]==0,abs(toy$MinorityClusterF1[2]-2/3)<1e-12,toy$MinorityClusterF1[3]==1,abs(toy$MinorityDetectionF1_all_extra_v2[3]-2/3)<1e-12)
 write(toy,'fixed_label_counterexamples.csv','evidence')
 # Prepare a factor-coverage diagnostic proposal, NOT new fitting.
 manifests<-list();proposal<-list();resource<-list()
 for(st in c('architecture','stress')) {
  m<-read.csv(file.path(project,'Scientific Reports Manuscript/github_results/simulation_design',paste0(st,'.csv')))
  m<-m[!duplicated(m$cell_id),]
  if(st=='architecture')m<-m[m$n==1000&m$p==5&m$separation==6&m$main_structure=='elliptical',]
  stopifnot(nrow(m)==if(st=='architecture')5 else 30)
  m$original_task_id<-m$task_id;m$task_id<-NULL;m$rep_start<-1;m$rep_end<-20
  m$proposal_id<-paste0('C3-',st,'-',m$cell_id);m$execution_status<-'NOT EXECUTED; specific author approval required'
  proposal[[st]]<-m
  a<-original[[st]];a<-a[a$cell_id%in%m$cell_id,]
  resource[[st]]<-do.call(rbind,lapply(m$cell_id,function(cid){v<-a$end_to_end_seconds[a$cell_id==cid];data.frame(study=st,cell_id=cid,n_historical_fits=length(v),historical_mean_fit_elapsed=mean(v),historical_median_fit_elapsed=median(v),planned_datasets=20,planned_method_fits=60,projected_elapsed_sum_seconds=mean(v)*60)}))
 }
 write(do.call(rbind,proposal),'targeted_validation_proposal_NOT_EXECUTED.csv');write(do.call(rbind,resource),'targeted_validation_resource_estimate.csv')
 # Scalar-only audit outputs remain compact; do not duplicate the full rows.
 saveRDS(list(summaries=sumtab,scenarios=scen,paired=pairtab,diagnostics=do.call(rbind,diagnostic),counts=do.call(rbind,counts)),file.path(c3,'evidence/compact_otrimle_audit.rds'))
 checks<-data.frame(check=c('scientific source hashes unchanged','Tier1 record and pairing integrity','original RDS versus submitted gzip scalar exports','original headline pooled means','all local auxiliary RDS have scalar-only schemas','fixed-label metric distinctions'),status='PASS',detail=c(paste(nrow(inv),'inputs'),'400 OTRIMLE records, 4 scenarios; exact shared data/starting seed checks; 1200 three-method records','Architecture, stress and negative-control rows, keys and available metrics agree','24 method/metric targets agree within 1e-11',paste(nrow(do.call(rbind,availability)),'objects inspected'),'No empirical simulation or fitting; full precision and extra-set handling verified'))
 write(checks,'validation_checks.csv','evidence');writeLines(capture.output(sessionInfo()),file.path(c3,'evidence/R_session_info.txt'))
 print(do.call(rbind,diagnostic));cat('COMMENT3_READONLY_AUDIT_OK\n')
}
main()

# Companion non-fitting audits refine the retention-only diagnostic table.
source(file.path(code,'audit_component_capacity.R'))
source(file.path(code,'audit_sufficient_rates.R'))
source(file.path(code,'test_origin_bounds.R'))
source(file.path(code,'audit_export_columns.R'))
# Keep the compact object explicitly current while retaining the coarser bounds.
z<-readRDS(file.path(c3,'evidence/compact_otrimle_audit.rds'))
z$retention_only_diagnostics<-z$diagnostics
z$diagnostics<-read.csv(file.path(c3,'tables/otrimle_macro_minority_diagnostics.csv'))
saveRDS(z,file.path(c3,'evidence/compact_otrimle_audit.rds'))
