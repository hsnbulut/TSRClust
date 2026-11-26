# Read-only summaries of frozen results. No fitting, simulated observations or jobs.
stopifnot(as.character(getRversion())=='4.5.2')
a<-grep('^--file=',commandArgs(),value=TRUE)[1]
C<-dirname(dirname(normalizePath(sub('^--file=','',a))));R<-dirname(dirname(C));P<-dirname(R)
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
inv<-read.csv(file.path(C,'evidence/source_inventory.csv'))
for(i in seq_len(nrow(inv)))stopifnot(sha(if(nzchar(inv$project_relative[i]))file.path(P,inv$project_relative[i]) else archive)==inv$sha256[i])
bind<-function(x)if(length(x))do.call(rbind,x) else data.frame()
wr<-function(x,n,w='tables')write.csv(x,file.path(C,w,n),row.names=FALSE,na='NA')
interval<-function(v,binary=FALSE){
 v<-v[is.finite(v)];n<-length(v)
 if(!n)return(c(n_finite=0,events=NA,mean=NA,mcse=NA,lower=NA,upper=NA))
 se<-if(n>1)sd(v)/sqrt(n) else NA_real_;ci<-if(n>1)mean(v)+c(-1,1)*qt(.975,n-1)*se else c(NA,NA)
 if(binary){stopifnot(all(v%in%0:1));ci<-binom.test(sum(v),n)$conf.int}
 c(n_finite=n,events=if(binary)sum(v) else NA,mean=mean(v),mcse=se,lower=ci[1],upper=ci[2])
}
pairci<-function(a,b,binary=FALSE){
 ok<-is.finite(a)&is.finite(b);a<-a[ok];b<-b[ok];n<-length(a);z<-interval(a-b);up<-down<-NA
 if(binary&&n){up<-sum(a==1&b==0);down<-sum(a==0&b==1);u<-binom.test(up,n,conf.level=.975)$conf.int;d<-binom.test(down,n,conf.level=.975)$conf.int;z[c('lower','upper')]<-c(u[1]-d[2],u[2]-d[1])}
 c(z,n_up=up,n_down=down)
}
main<-function(){
 t<-tempfile('r2c2_archive_');dir.create(t);on.exit(unlink(t,recursive=TRUE))
 untar(archive,files='replications.csv',exdir=t);all<-read.csv(file.path(t,'replications.csv'))
 design<-read.csv(file.path(P,'revision_scientific_reports/config/design.csv'));design<-design[design$tier==1,]
 tasks<-read.csv(file.path(P,'revision_scientific_reports/config/tasks_tier1.csv'))
 stopifnot(nrow(all)==98800,all(all$status=='OK'),length(unique(all$cell_id))==282,nrow(tasks)==9880,!anyDuplicated(all[c('cell_id','replication')]))
 keys<-paste(all$cell_id,all$replication);tk<-unlist(lapply(seq_len(nrow(tasks)),function(i)paste(tasks$cell_id[i],seq.int(tasks$rep_start[i],tasks$rep_end[i]))));stopifnot(!anyDuplicated(tk),setequal(keys,tk))
 idx<-match(all$cell_id,design$cell_id)
 for(nm in intersect(names(all),names(design)))stopifnot(identical(is.na(all[[nm]]),is.na(design[[nm]][idx])),all(as.character(all[[nm]])==as.character(design[[nm]][idx]),na.rm=TRUE))
 env<-new.env();sys.source(file.path(P,'revision_scientific_reports/R/common.R'),env)
 for(s in c('data','algorithm','null')){
  stream<-match(s,c('data','algorithm','null'))-1L;rep<-if(s=='data')rep(1L,nrow(all)) else all[[paste0(s,'_repeat')]]
  stopifnot(all(all[[paste0(s,'_seed')]]==mapply(env$seed_id,all$dataset_id,all$replication,rep,MoreArgs=list(stream=stream))))
 }
 x<-all[all$method=='TSRClust',];stopifnot(nrow(x)==94800,all(x$structure_alpha==.05),all(x$any_native_micro==as.integer(x$n_accepted_micro>0)))
 baselineC<-x$experiment=='C'&x$density==2&x$novelty==.995&x$minpts_multiplier==2&x$B%in%c(49,99,199)
 z<-x[x$experiment%in%c('A-END','A-POWER','A-COND')|baselineC,]
 stopifnot(nrow(z)==64000,length(unique(z$cell_id))==80)
 meta<-c('experiment','scenario','cell_id','dataset_id','rule','reference','mode','n','p','q','minority_ratio','noise_type','requested_alpha','nR','B','structure_alpha','density','novelty','minpts_multiplier')
 metrics<-c('any_native_micro','AnyMicroDiscovery','MinorityClusterF1','MinorityDetectionF1_all_extra_v2','SpuriousExtraClusters','n_accepted_micro','residual_assignment_rate','end_to_end_seconds','micro_seconds','macro_seconds','global_p','FalseExtraAssignmentRate')
 binary<-c('any_native_micro','AnyMicroDiscovery')
 cells<-reasons<-list()
 for(cid in unique(z$cell_id)){
  a<-z[z$cell_id==cid,];m<-a[1,meta];rownames(m)<-NULL
  for(v in metrics)cells[[length(cells)+1]]<-cbind(m,n_records=nrow(a),metric=v,as.data.frame(as.list(interval(a[[v]],v%in%binary))),median=if(any(is.finite(a[[v]])))median(a[[v]],na.rm=TRUE) else NA_real_,uncertainty=if(v%in%binary)'Pointwise exact binomial95% within one fixed cell' else 'Pointwise t Monte Carlo95%;timing includes execution/resource variation')
  for(reason in unique(a$rejection_reason))reasons[[length(reasons)+1]]<-cbind(m,reason=reason,n_records=nrow(a),count=sum(a$rejection_reason==reason),rate=mean(a$rejection_reason==reason))
 }
 wr(bind(cells),'cell_summaries.csv');wr(bind(reasons),'rejection_reasons.csv')
 availability<-all[!duplicated(all$cell_id),names(design)];availability$available_records<-vapply(availability$cell_id,function(cid)sum(all$cell_id==cid),integer(1));wr(availability,'experiment_inventory.csv')
 proof<-pairs<-conditional<-breps<-changes<-list();counter<-0L
 fields<-c('n','p','main_k','separation','main_structure','main_distribution','main_weights','residual_total','minority_ratio','q','minority_shape','minority_sd','minority_proximity','noise_type','mode','nR','data_seed','algorithm_seed','null_seed','algorithm_repeat','null_repeat','algorithm_version','metric_version','reference','novelty','density','minpts_multiplier','B','nstart','regime','kmax','structure_alpha','restr_fact','eps_grid_size','requested_alpha')
 compare<-function(a,b,family,vary='rule'){
  a<-a[order(a$replication),];b<-b[order(b$replication),];stopifnot(nrow(a)==nrow(b),identical(a$replication,b$replication),identical(a$dataset_id,b$dataset_id))
  fs<-setdiff(fields,vary);stopifnot(all(vapply(fs,function(v)identical(a[[v]],b[[v]]),logical(1))))
  stopifnot(identical(a$n_residual,b$n_residual),identical(a$macro_objective,b$macro_objective),identical(a$macro_code,b$macro_code))
  if(vary=='rule'&&a$rule[1]=='pvalue_only'&&b$rule[1]=='legacy_rule')stopifnot(identical(a$global_p,b$global_p))
  counter<<-counter+1L;id<-sprintf('R2C2-%02d',counter)
  row<-data.frame(contrast_id=id,family=family,scenario=a$scenario[1],q=a$q[1],p=a$p[1],minority_ratio=a$minority_ratio[1],noise_type=a$noise_type[1],requested_alpha=a$requested_alpha[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],new_rule=a$rule[1],baseline_rule=b$rule[1],new_B=a$B[1],baseline_B=b$B[1])
  proof[[length(proof)+1]]<<-cbind(row,n_pairs=nrow(a),dataset_id=a$dataset_id[1],all_noncontrast_fields_equal=TRUE,all_starting_seeds_equal=TRUE,macro_objective_code_residual_equal=TRUE,data_hashes_exported=FALSE)
  out<-list()
  for(v in metrics)out[[length(out)+1]]<-cbind(row,metric=v,mean_new=if(any(is.finite(a[[v]])))mean(a[[v]],na.rm=TRUE) else NA_real_,mean_baseline=if(any(is.finite(b[[v]])))mean(b[[v]],na.rm=TRUE) else NA_real_,as.data.frame(as.list(pairci(a[[v]],b[[v]],v%in%binary))),direction='New minus baseline',uncertainty=if(v%in%binary)'Conservative exact paired95% bound from two97.5% CP discordance margins;pointwise' else 'Paired t Monte Carlo95%;timing not a randomized execution comparison')
  ch<-cbind(row,n_pairs=nrow(a),native_up=sum(a$any_native_micro==1&b$any_native_micro==0),native_down=sum(a$any_native_micro==0&b$any_native_micro==1),both_native=sum(a$any_native_micro==1&b$any_native_micro==1),accepted_count_changed=sum(a$n_accepted_micro!=b$n_accepted_micro),residual_fraction_changed=sum(abs(a$residual_assignment_rate-b$residual_assignment_rate)>1e-12),selected_eps_changed=sum((is.na(a$selected_eps)!=is.na(b$selected_eps))|(!is.na(a$selected_eps)&!is.na(b$selected_eps)&a$selected_eps!=b$selected_eps)),cluster_F1_changed=sum(abs(a$MinorityClusterF1-b$MinorityClusterF1)>1e-12,na.rm=TRUE),detection_F1_changed=sum(abs(a$MinorityDetectionF1_all_extra_v2-b$MinorityDetectionF1_all_extra_v2)>1e-12,na.rm=TRUE))
  if(vary=='rule'&&a$rule[1]=='pvalue_only'&&b$rule[1]=='legacy_rule'){
   up<-a$any_native_micro==1&b$any_native_micro==0;stopifnot(all(a$global_p[up]==.04),all(b$rejection_reason[up]=='candidate_max'),!any(a$any_native_micro<b$any_native_micro))
   ch$native_up_global_p_004<-sum(up&a$global_p==.04);ch$global_p_equal<-TRUE
  }
  changes[[length(changes)+1]]<<-ch
  bind(out)
 }
 for(ex in c('A-END','A-POWER'))for(sc in unique(z$scenario[z$experiment==ex])){
  a<-z[z$experiment==ex&z$scenario==sc,];pairs[[length(pairs)+1]]<-compare(a[a$rule=='pvalue_only',],a[a$rule=='legacy_rule',],paste(ex,'matched rules'))
 }
 for(sc in unique(z$scenario[z$experiment=='A-COND'])){
  a<-z[z$experiment=='A-COND'&z$scenario==sc,]
  for(rs in list(c('pvalue_only','legacy_rule'),c('recalibrated_rule','legacy_rule'),c('recalibrated_rule','pvalue_only')))conditional[[length(conditional)+1]]<-compare(a[a$rule==rs[1],],a[a$rule==rs[2],],'A-COND distinct uniform residual-only')
 }
 for(sc in unique(z$scenario[z$experiment=='C'])){
  a<-z[z$experiment=='C'&z$scenario==sc,]
  for(b in c(99,199))breps[[length(breps)+1]]<-compare(a[a$B==b,],a[a$B==49,],'C matched B sensitivity','B')
 }
 wr(bind(pairs),'paired_legacy_vs_pvalue.csv');wr(bind(conditional),'paired_conditional_rules.csv');wr(bind(breps),'paired_B_vs49.csv');wr(bind(proof),'pairing_proofs.csv','evidence');wr(env$bind_rows(changes),'decision_change_counts.csv')
 # Reconcile approved C1 cell estimates and all overlapping paired contrasts.
 old<-read.csv(file.path(R,'reviewer2/comment1/tables/false_micro_discovery_by_noise.csv'))
 for(i in seq_len(nrow(old))){a<-x[x$cell_id==old$cell_id[i],];stopifnot(sum(a$any_native_micro)==old$native_events[i],sum(a$AnyMicroDiscovery)==old$external_events[i])}
 oc<-read.csv(file.path(R,'reviewer2/comment1/tables/paired_rule_differences.csv'))
 new<-rbind(bind(pairs),bind(conditional));nc<-0L
 for(i in seq_len(nrow(oc))){v<-oc$metric[i];if(!v%in%metrics)next
  a<-new[new$new_cell==oc$new_cell[i]&new$baseline_cell==oc$baseline_cell[i]&new$metric==v,];stopifnot(nrow(a)==1,abs(a$mean-oc$mean[i])<1e-12);nc<-nc+1L
 }
 accepted<-x$any_native_micro==1
 stopifnot(all(abs(x$global_p[accepted&x$rule=='legacy_rule']-1/(x$B[accepted&x$rule=='legacy_rule']+1))<1e-14))
 wr(data.frame(check=c('Immutable source inputs','Complete frozen Tier1 records/configurations/tasks','Design fields and all deterministic seeds','80 targeted cells and12 metrics','48 paired groups with noncontrast and Stage-I identity','Approved C1 event frequencies and overlapping paired differences','Legacy accepted global p equals minimum attainable p','No new fits or generated observations'),status='PASS',detail=c('235 SHA256 unchanged','98800 OK;282 configurations;9880 blocks','Every common design field andseed formula matches','64000 records;960 summaries','28 A-END/A-POWER;12 A-COND;8 B contrasts',paste('81 cells;',nc,'overlapping differences exact'),'All94800 TSR records checked','Read-only aggregation and boundary fixtures only')),'analysis_checks.csv','evidence')
 saveRDS(z,file.path(C,'evidence/audited_scalar_records.rds'))
 writeLines(capture.output(sessionInfo()),file.path(C,'evidence/R_session_info.txt'))
 cat('R2C2_ANALYSIS_OK:80 cells;48paired groups;960cell summaries;C1 unchanged\n')
}
main()
