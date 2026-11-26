# Aggregate count identifiability from PER-FIT P/R and recorded rate identities.
# Never invert ARI/F1 or reconstruct observation-level assignments. No fitting.
f<-grep('^--file=',commandArgs(),value=TRUE)[1];c3<-dirname(dirname(normalizePath(sub('^--file=','',f))));project<-dirname(dirname(dirname(c3)))
source(file.path(c3,"code/origin_bounds.R"))
run_sufficient_audit<-function(){
 work<-tempfile('c3_sufficient_');dir.create(work);on.exit(unlink(work,recursive=TRUE))
 chunks<-list()
 for(st in c('architecture','stress','negative_control')){
  x<-readRDS(file.path(project,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined',paste0(st,'_replications.rds')));x<-x[x$method=='OTRIMLE_fixedK',]
  x$source_family<-'Original submitted';x$cohort<-if(st=='architecture')ifelse(x$q>0,'Original architecture M1-M4','Original architecture N0') else paste('Original',st)
  x$pred_zero_count_for_audit<-x$estimated_noise_proportion*x$n
  chunks[[st]]<-x
 }
 untar(Sys.getenv('TSR_TIER1_ARCHIVE'),files='replications.csv',exdir=work)
 x<-read.csv(file.path(work,'replications.csv'));x<-x[x$method=='OTRIMLE_fixedK',];x$source_family<-'Tier1 revision';x$cohort<-paste0('Tier1 ',x$scenario);x$pred_zero_count_for_audit<-x$n_residual;chunks[['tier1']]<-x
 numeric_check<-list();records<-list();j<-0
 nearest_integer<-function(v,name,nmax){
  stopifnot(all(is.finite(v)),all(v>= -1e-7),all(v<=nmax+1e-7))
  error<-max(abs(v-round(v)));stopifnot(error<1e-7)
  numeric_check[[length(numeric_check)+1]]<<-data.frame(quantity=name,n=length(v),max_integer_residual=error,tolerance=1e-7,status='PASS')
  as.integer(round(v))
 }
 for(source in names(chunks)){
  a<-chunks[[source]];n<-a$n
  # These population cardinalities follow the saved generator's exact lines:
  # main_total = n-round(n*residual_total); minority_total=round(n*minority_ratio).
  M<-n-round(n*a$residual_total);U<-round(n*a$minority_ratio);N<-n-M-U
  stopifnot(all(M>0),all(U>=0),all(N>=0),all(U[a$q==0]==0))
  Vcount<-a$QHat-a$SpuriousExtraClusters
  TP<-nearest_integer(ifelse(U>0,a$MinorityRecall*U,0),paste(source,'TP from recall'),n)
  positive<-Vcount>0
  stopifnot(all((TP>0)==positive),all(a$MinorityPrecision[positive]>0))
  predV<-numeric(nrow(a));predV[positive]<-TP[positive]/a$MinorityPrecision[positive]
  predV<-nearest_integer(predV,paste(source,'V union size from per-fit precision'),n)
  FP<-predV-TP;FN<-U-TP;stopifnot(all(FP>=0),all(FN>=0),all(predV[!positive]==0))
  # Cross-check against F1, without using F1 to infer any count.
  defined<-U>0
  f1check<-2*TP[defined]/(2*TP[defined]+FP[defined]+FN[defined])
  stopifnot(all(abs(f1check-a$MinorityDetectionF1[defined])<1e-11))
  macroP<-nearest_integer(a$MajorRetention*M,paste(source,'macro P'),n);nonret<-M-macroP
  FPE<-nearest_integer(a$FalseExtraAssignmentRate*(M+N),paste(source,'nonminority E'),n)
  noise0<-nearest_integer(ifelse(N>0,a$NoiseRecall*N,0),paste(source,'noise zero'),n)
  all0<-nearest_integer(a$pred_zero_count_for_audit,paste(source,'all predicted zero'),n)
  stopifnot(all(FPE>=FP),all(noise0<=N),all(all0>=noise0),all(FN>=0))
  # x_E is the true macro count in E. Known identities give valid integer bounds.
  # macro0=nonret-x_E; minority0=all0-noise0-nonret+x_E.
  # Every unmatched extra cluster needs >=1 member; any nonminority count
  # shortfall must be made up by true minorities outside the V union.
  bounds<-origin_bounds(M,U,N,TP,FP,FN,macroP,FPE,noise0,all0,a$SpuriousExtraClusters)
  LE<-bounds$LE;UE<-bounds$UE
  if(source=='tier1'){
   direct<-nearest_integer(a$MacroFalsePositiveReassignment*M,'Tier1 direct macro E',n)
   stopifnot(all(direct>=LE),all(direct<=UE));LE<-UE<-direct
  }
  LV<-pmax(0,FP-(N-noise0),LE-(FPE-FP));UV<-pmin(FP,UE)
  # V subset E, with FPE-FP nonminority observations in E minus V.
  LE<-pmax(LE,LV);UE<-pmin(UE,UV+(FPE-FP))
  LV<-pmax(LV,LE-(FPE-FP));UV<-pmin(UV,UE)
  stopifnot(all(LE<=UE),all(LV<=UV),all(LV>=0),all(UV<=M),all(UE<=M))
  empty<-Vcount==0;stopifnot(all(LV[empty]==0),all(UV[empty]==0))
  j<-j+1;records[[j]]<-data.frame(source_family=a$source_family,cohort=a$cohort,cell_id=a$cell_id,replication=a$replication,n=n,n_macro=M,n_minority=U,n_noise=N,TP_V=TP,FP_V=FP,FN_V=FN,macro_P=macroP,noise_0=noise0,predicted_0=all0,FP_E=FPE,macro_V_lower_count=LV,macro_V_upper_count=UV,macro_E_lower_count=LE,macro_E_upper_count=UE,macro_V_lower_rate=LV/M,macro_V_upper_rate=UV/M,macro_E_lower_rate=LE/M,macro_E_upper_rate=UE/M)
 }
 all<-do.call(rbind,records);saveRDS(all,file.path(c3,'evidence/identified_aggregate_counts.rds'),compress='xz')
 write.csv(do.call(rbind,numeric_check),file.path(c3,'evidence/sufficient_rate_integer_checks.csv'),row.names=FALSE)
 results<-list();counts<-list();j<-0
 for(cohort in unique(all$cohort)){
  a<-all[all$cohort==cohort,]
  counts[[length(counts)+1]]<-data.frame(source_family=a$source_family[1],cohort=cohort,n_records=nrow(a),n_TP_FP_FN_identified=nrow(a),mean_TP_V=mean(a$TP_V),mean_FP_V=mean(a$FP_V),mean_FN_V=mean(a$FN_V),TP_V_sum=sum(a$TP_V),FP_V_sum=sum(a$FP_V),FN_V_sum=sum(a$FN_V),basis='Aggregate counts deduced from per-fit precision/recall and deterministic generator cardinalities; integer residual checks <1e-7; NO ARI/F1 inversion; NO individual memberships reconstructed')
  for(target in c('Minority','AnyExtra')){
   label<-if(target=='Minority')'V' else 'E';lower<-a[[paste0('macro_',label,'_lower_rate')]];upper<-a[[paste0('macro_',label,'_upper_rate')]];exact<-lower==upper
   j<-j+1;results[[j]]<-data.frame(source_family=a$source_family[1],cohort=cohort,target=paste0('MacroTo',target),denominator='All true macro observations per fit; cohort summary averages fit-specific fractions',n_records=nrow(a),n_exactly_identified=sum(exact),n_NOT_IDENTIFIABLE=sum(!exact),mean_lower_bound=mean(lower),mean_upper_bound=mean(upper),sum_lower_count=sum(a[[paste0('macro_',label,'_lower_count')]]),sum_upper_count=sum(a[[paste0('macro_',label,'_upper_count')]]),n_exact_zero=sum(exact&upper==0),n_definitely_positive=sum(lower>0),exact_point_rate_if_all_identified=if(all(exact))mean(upper) else NA_real_,status=if(all(exact))'IDENTIFIED AGGREGATE COUNTS; NOT SAVED MEMBERSHIPS' else 'NOT IDENTIFIABLE FOR ALL FITS: sufficient-rate bounds',basis='Per-fit P/R, retention, false-extra rate, noise recall, hard-zero total, deterministic cardinalities and E/V counts. All integer-compatible. No ARI/F1 inversion.')
  }
 }
 old<-file.path(c3,'tables/otrimle_macro_minority_diagnostics.csv')
 if(file.exists(old)&&!file.exists(file.path(c3,'tables/retention_only_promotion_bounds.csv')))file.copy(old,file.path(c3,'tables/retention_only_promotion_bounds.csv'))
 write.csv(do.call(rbind,results),old,row.names=FALSE);write.csv(do.call(rbind,counts),file.path(c3,'tables/identified_detection_totals.csv'),row.names=FALSE)
 amb<-all[all$macro_V_lower_count!=all$macro_V_upper_count|all$macro_E_lower_count!=all$macro_E_upper_count,]
 write.csv(amb,file.path(c3,'evidence/ambiguous_origin_count_ranges.csv'),row.names=FALSE)
 cat('SUFFICIENT_RATE_AGGREGATE_AUDIT_OK: no ARI/F1 inversion, no fitting; rows',nrow(all),'\n');print(do.call(rbind,results));print(do.call(rbind,counts))
}
run_sufficient_audit()
