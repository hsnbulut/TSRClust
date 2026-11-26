# Read-only Reviewer2 Comment1 audit. No generated observations, fits or jobs.
stopifnot(as.character(getRversion())=='4.5.2')
arg<-grep('^--file=',commandArgs(),value=TRUE)[1];C<-dirname(dirname(normalizePath(sub('^--file=','',arg))));R<-dirname(dirname(C));P<-dirname(R)
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
inv<-read.csv(file.path(C,'evidence/source_inventory.csv'))
for(i in seq_len(nrow(inv))){p<-if(nzchar(inv$project_relative[i]))file.path(P,inv$project_relative[i]) else archive;stopifnot(sha(p)==inv$sha256[i])}
bind<-function(x)if(length(x))do.call(rbind,x) else data.frame()
wr<-function(x,name,where='tables')write.csv(x,file.path(C,where,name),row.names=FALSE,na='NA')
summary_ci<-function(v,binary=FALSE){
 v<-v[is.finite(v)];n<-length(v)
 if(!n)return(c(n_finite=0,events=NA,mean=NA,mcse=NA,lower=NA,upper=NA))
 se<-if(n>1)sd(v)/sqrt(n) else NA_real_;ci<-if(n>1)mean(v)+c(-1,1)*qt(.975,n-1)*se else c(NA,NA)
 if(binary){stopifnot(all(v%in%c(0,1)));ci<-binom.test(sum(v),n)$conf.int}
 c(n_finite=n,events=if(binary)sum(v) else NA,mean=mean(v),mcse=se,lower=ci[1],upper=ci[2])
}
paired_ci<-function(a,b,binary=FALSE){
 ok<-is.finite(a)&is.finite(b);a<-a[ok];b<-b[ok];n<-length(a);z<-summary_ci(a-b);up<-down<-NA
 if(binary&&n){up<-sum(a==1&b==0);down<-sum(a==0&b==1);u<-binom.test(up,n,conf.level=.975)$conf.int;d<-binom.test(down,n,conf.level=.975)$conf.int;z[c('lower','upper')]<-c(u[1]-d[2],u[2]-d[1])}
 c(z,n_up=up,n_down=down)
}
work<-tempfile('r2c1_');dir.create(work)
main<-function(){
 on.exit(unlink(work,recursive=TRUE))
 untar(archive,files=c('replications.csv','micro_counts.csv'),exdir=work);all<-read.csv(file.path(work,'replications.csv'))
 design<-read.csv(file.path(P,'revision_scientific_reports/config/design.csv'));d1<-design[design$tier==1,]
 tasks<-read.csv(file.path(P,'revision_scientific_reports/config/tasks_tier1.csv'))
 stopifnot(nrow(all)==98800,all(all$status=='OK'),length(unique(all$cell_id))==282,nrow(d1)==282,nrow(tasks)==9880,!anyDuplicated(all[c('cell_id','replication')]))
 stopifnot(all(vapply(split(all,all$cell_id),function(a)nrow(a)==design$reps[match(a$cell_id[1],design$cell_id)],logical(1))))
 key<-paste(all$cell_id,all$replication);taskkey<-unlist(lapply(seq_len(nrow(tasks)),function(i)paste(tasks$cell_id[i],seq.int(tasks$rep_start[i],tasks$rep_end[i]))));stopifnot(!anyDuplicated(taskkey),setequal(key,taskkey))
 idx<-match(all$cell_id,design$cell_id)
 for(nm in intersect(names(design),names(all))){a<-all[[nm]];b<-design[[nm]][idx];stopifnot(identical(is.na(a),is.na(b)),all(as.character(a)==as.character(b),na.rm=TRUE))}
 env<-new.env();sys.source(file.path(P,'revision_scientific_reports/R/common.R'),env)
 stopifnot(all(all$data_seed==mapply(env$seed_id,all$dataset_id,all$replication,MoreArgs=list(stream=0L))),all(all$algorithm_seed==mapply(env$seed_id,all$dataset_id,all$replication,all$algorithm_repeat,MoreArgs=list(stream=1L))),all(all$null_seed==mapply(env$seed_id,all$dataset_id,all$replication,all$null_repeat,MoreArgs=list(stream=2L))))
 x<-all[all$q==0 & all$method=='TSRClust',];stopifnot(nrow(x)==63800,all(is.na(x$MinorityClusterF1)),all(x$QHat==x$SpuriousExtraClusters),all(x$AnyMicroDiscovery==as.integer(x$QHat>0)),all(x$any_native_micro==as.integer(x$n_accepted_micro>0)))
 # requested_alpha is the TCLUST trimming input,not the validation level.
 stopifnot(all(x$structure_alpha==.05),all(x$minpts_multiplier>=1))
 x$AnySpuriousExtra<-as.integer(x$SpuriousExtraClusters>0)
 x$accepted_assignment_count<-round(x$residual_assignment_rate*x$n_residual)
 stopifnot(max(abs(x$residual_assignment_rate*x$n_residual-x$accepted_assignment_count))<1e-7)
 x$accepted_assignment_fraction_all<-x$accepted_assignment_count/ifelse(x$mode=='residual_only',x$n_residual,x$n)
 x$derived_support<-pmax(5,x$minpts_multiplier*x$p)
 x$generated_macro_count<-ifelse(x$mode=='residual_only',0,x$n-round(x$n*x$residual_total))
 x$generated_noise_count<-ifelse(x$mode=='residual_only',x$n_residual,x$n-x$generated_macro_count)
 # Without any true minority,all E assignments are false positives relative to truth.
 x$external_extra_count<-round(x$FalseExtraAssignmentRate*ifelse(x$mode=='residual_only',x$n_residual,x$n))
 stopifnot(max(abs(x$FalseExtraAssignmentRate*ifelse(x$mode=='residual_only',x$n_residual,x$n)-x$external_extra_count))<1e-7)
 x$macro_extra_count<-round(ifelse(x$generated_macro_count>0,x$MacroFalsePositiveReassignment*x$generated_macro_count,0))
 x$noise_extra_count<-x$external_extra_count-x$macro_extra_count
 stopifnot(all(x$noise_extra_count>=0),all(x$noise_extra_count<=x$generated_noise_count))
 x$noise_to_external_extra_rate<-x$noise_extra_count/x$generated_noise_count
 metrics<-c('AnyMicroDiscovery','any_native_micro','SpuriousExtraClusters','n_accepted_micro','residual_assignment_rate','accepted_assignment_count','accepted_assignment_fraction_all','FalseExtraAssignmentRate','MacroFalsePositiveReassignment','noise_to_external_extra_rate','MajorRetention','n_residual')
 binary<-c('AnyMicroDiscovery','any_native_micro','AnySpuriousExtra')
 configs<-x[!duplicated(x$cell_id),intersect(names(design),names(x))];configs$n_records<-vapply(configs$cell_id,function(cid)sum(x$cell_id==cid),integer(1));wr(configs,'negative_configuration_inventory.csv')
 summaries<-list();overview<-list();reasons<-list()
 keys<-c('experiment','scenario','cell_id','dataset_id','rule','reference','mode','noise_type','n','p','nR','residual_total','requested_alpha','structure_alpha','B','density','novelty','minpts_multiplier','derived_support','regime','kmax')
 for(cid in unique(x$cell_id)){
  a<-x[x$cell_id==cid,];meta<-a[1,keys];row.names(meta)<-NULL;meta$n_records<-nrow(a)
  stopifnot(nrow(a)==design$reps[match(cid,design$cell_id)])
  for(m in metrics)summaries[[length(summaries)+1]]<-cbind(meta,metric=m,as.data.frame(as.list(summary_ci(a[[m]],m%in%binary))),uncertainty=if(m%in%binary)'Pointwise exact Clopper-Pearson95% for one fixed generating cell' else 'Pointwise t Monte Carlo95%;zero observed variance is not a population guarantee',source='Immutable Tier1 replications.csv')
  ev<-summary_ci(a$AnyMicroDiscovery,TRUE);nv<-summary_ci(a$any_native_micro,TRUE)
  overview[[length(overview)+1]]<-cbind(meta,external_events=ev['events'],external_rate=ev['mean'],external_lower=ev['lower'],external_upper=ev['upper'],native_events=nv['events'],native_rate=nv['mean'],native_lower=nv['lower'],native_upper=nv['upper'],event_disagreements=sum(a$AnyMicroDiscovery!=a$any_native_micro),mean_spurious=mean(a$SpuriousExtraClusters),mean_accepted_clusters=mean(a$n_accepted_micro),mean_residual_assigned=mean(a$residual_assignment_rate),mean_false_extra_assignment=mean(a$FalseExtraAssignmentRate),mean_macro_extra=if(any(is.finite(a$MacroFalsePositiveReassignment)))mean(a$MacroFalsePositiveReassignment,na.rm=TRUE) else NA_real_,mean_noise_extra=mean(a$noise_to_external_extra_rate),denominator='Frequency uses independent datasets within one fixed cell;residual fraction uses that fit\'s actual residual count;false-extra fraction uses all nonminority observations',source='Immutable Tier1 replications.csv')
  for(reason in unique(a$rejection_reason))reasons[[length(reasons)+1]]<-cbind(meta,reason=reason,count=sum(a$rejection_reason==reason),frequency=mean(a$rejection_reason==reason))
 }
 wr(bind(summaries),'negative_cell_metric_summaries.csv');wr(bind(overview),'false_micro_discovery_by_noise.csv');wr(bind(reasons),'negative_rejection_reasons.csv')
 wr(x[,c(keys,'replication','data_seed','algorithm_seed','null_seed','AnyMicroDiscovery','any_native_micro','n_accepted_micro','n_residual','accepted_assignment_count','external_extra_count','macro_extra_count','noise_extra_count')],'derived_negative_assignment_counts.csv')
 # Genuine pair proofs:every DGP field and all three starting seed streams.
 dgp<-c('n','p','main_k','separation','main_structure','main_distribution','main_weights','residual_total','minority_ratio','q','minority_shape','minority_sd','minority_proximity','noise_type','mode','nR')
 streams<-c('data_seed','algorithm_seed','null_seed','algorithm_repeat','null_repeat','algorithm_version','metric_version')
 settings<-c('reference','novelty','density','minpts_multiplier','B','nstart','regime','kmax','structure_alpha','restr_fact','eps_grid_size','requested_alpha')
 pairs<-proof<-list();counter<-0L
 compare<-function(a,b,family,whichmetrics=metrics){
  a<-a[order(a$replication),];b<-b[order(b$replication),];stopifnot(nrow(a)==nrow(b),identical(a$replication,b$replication),all(a$dataset_id==b$dataset_id))
  stopifnot(all(vapply(c(dgp,streams,settings),function(nm)identical(a[[nm]],b[[nm]]),logical(1))))
  stopifnot(identical(a$n_residual,b$n_residual),identical(a$macro_objective,b$macro_objective),identical(a$macro_code,b$macro_code))
  counter<<-counter+1L;id<-sprintf('R2C1-%02d',counter)
  proof[[length(proof)+1]]<<-data.frame(contrast_id=id,family=family,scenario=a$scenario[1],new_rule=a$rule[1],baseline_rule=b$rule[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],n_pairs=nrow(a),dataset_id=a$dataset_id[1],all_DGP_settings_equal=TRUE,all_nonrule_algorithm_settings_equal=TRUE,all_starting_seeds_equal=TRUE,macro_objective_code_and_residual_equal=TRUE,interpretation='Exact pairing under frozen generator/seed contract;actual X hashes and realized null draws not exported. Recalibrated changes both null generator and own-grid recomputation.')
  for(m in whichmetrics)pairs[[length(pairs)+1]]<<-cbind(data.frame(contrast_id=id,family=family,scenario=a$scenario[1],noise_type=a$noise_type[1],p=a$p[1],requested_alpha=a$requested_alpha[1],nR=a$nR[1],new_rule=a$rule[1],baseline_rule=b$rule[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],metric=m,mean_new=if(any(is.finite(a[[m]])))mean(a[[m]],na.rm=TRUE) else NA_real_,mean_baseline=if(any(is.finite(b[[m]])))mean(b[[m]],na.rm=TRUE) else NA_real_),as.data.frame(as.list(paired_ci(a[[m]],b[[m]],m%in%binary))),direction='New minus baseline',uncertainty=if(m%in%binary)'Conservative exact paired95% bound from two97.5% Clopper-Pearson discordance margins;pointwise' else 'Pointwise paired t Monte Carlo95%;not multiplicity-adjusted',source='Immutable Tier1 replications.csv')
 }
 for(sc in unique(x$scenario[x$experiment=='A-END'])){a<-x[x$experiment=='A-END'&x$scenario==sc,];compare(a[a$rule=='pvalue_only',],a[a$rule=='legacy_rule',],'A-END matched rules')}
 for(sc in unique(x$scenario[x$experiment=='A-COND'])){a<-x[x$experiment=='A-COND'&x$scenario==sc,];compare(a[a$rule=='pvalue_only',],a[a$rule=='legacy_rule',],'A-COND matched rules');compare(a[a$rule=='recalibrated_rule',],a[a$rule=='legacy_rule',],'A-COND matched rules');compare(a[a$rule=='recalibrated_rule',],a[a$rule=='pvalue_only',],'A-COND matched rules')}
 wr(bind(pairs),'paired_rule_differences.csv');wr(bind(proof),'decision_pairing_proofs.csv','evidence')
 # Original submitted negative controls,never pooled with revision populations.
 path<-file.path(P,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined/negative_control_replications.rds');orig<-readRDS(path)
 ex<-read.csv(gzfile(file.path(P,'Scientific Reports Manuscript/github_results/raw_results/negative_control_replications.csv.gz')))
 stopifnot(nrow(orig)==80000,all(orig$status=='OK'),identical(names(orig),names(ex)))
 for(nm in names(orig)){
  a<-orig[[nm]];b<-ex[[nm]]
  if(is.character(a)){a[a=='']<-NA_character_;b<-as.character(b);b[b=='']<-NA_character_}
  stopifnot(identical(is.na(a),is.na(b)))
  if(is.numeric(a))stopifnot(all(abs(a-b)<1e-9,na.rm=TRUE)) else stopifnot(all(as.character(a)==as.character(b),na.rm=TRUE))
 }
 o<-orig[orig$method=='TSRClust_validated',];u<-orig[orig$method=='TSRClust_unvalidated',];stopifnot(nrow(o)==10000,nrow(u)==10000,all(o$q==0),all(o$null_reps==49),all(o$requested_k==2),all(o$nstart==500),all(o$requested_c==20))
 # With the strict maximum gate,positive score and p=1/(B+1) identify any native acceptance.
 o$native_event_derived<-as.integer(o$micro_score>0 & abs(o$micro_p_value-1/(o$null_reps+1))<1e-12)
 stopifnot(!anyNA(o$native_event_derived),sum(o$AnyMicroDiscovery)==3,sum(o$native_event_derived)==3,all(o$native_event_derived==o$AnyMicroDiscovery),sum(u$AnyMicroDiscovery)==10000)
 ocells<-o[!duplicated(o$cell_id),c('cell_id','n','p','main_k','separation','main_structure','main_distribution','main_weights','residual_total','q','minority_ratio','noise_type','requested_k','requested_alpha','requested_c','nstart','null_reps')]
 os<-opool<-omatch<-list()
 for(method in unique(orig$method)){
  a<-orig[orig$method==method,];e<-summary_ci(a$AnyMicroDiscovery,TRUE)
  opool[[length(opool)+1]]<-data.frame(method=method,n_records=nrow(a),n_cells=length(unique(a$cell_id)),events=e['events'],frequency=e['mean'],binomial95_lower=e['lower'],binomial95_upper=e['upper'],mean_spurious=mean(a$SpuriousExtraClusters),source=path,interval_scope='Reproduced conventional pooled binomial interval;heterogeneous fixed cells are not established iid Bernoulli trials. Cell-specific intervals below retain fixed-design meaning.')
 }
 for(cid in unique(o$cell_id)){
  a<-o[o$cell_id==cid,];b<-u[u$cell_id==cid,];a<-a[order(a$replication),];b<-b[order(b$replication),]
  stopifnot(nrow(a)==1000,identical(a$replication,b$replication),identical(a$seed,b$seed),all(vapply(dgp[!dgp%in%c('mode','nR')],function(nm)identical(a[[nm]],b[[nm]]),logical(1))))
  for(method in c('TSRClust_validated','TSRClust_unvalidated')){
   z<-if(method=='TSRClust_validated')a else b;ev<-summary_ci(z$AnyMicroDiscovery,TRUE)
   os[[length(os)+1]]<-cbind(ocells[ocells$cell_id==cid,],method=method,n_records=nrow(z),events=ev['events'],frequency=ev['mean'],lower=ev['lower'],upper=ev['upper'],mean_spurious=mean(z$SpuriousExtraClusters),native_event_derived=if(method=='TSRClust_validated')sum(z$native_event_derived) else NA,accepted_count_observed=NA_real_,assignment_fraction_observed=NA_real_,source=path)
  }
  zz<-paired_ci(a$AnyMicroDiscovery,b$AnyMicroDiscovery,TRUE);omatch[[length(omatch)+1]]<-cbind(data.frame(cell_id=cid,noise_type=a$noise_type[1],p=a$p[1],requested_alpha=a$requested_alpha[1],mean_new=mean(a$AnyMicroDiscovery),mean_baseline=mean(b$AnyMicroDiscovery)),as.data.frame(as.list(zz)),direction='Validated minus unvalidated',interpretation='Matched dataset and common Stage-I fit;not a pure Monte Carlo ablation:unvalidated uses fixed90th-percentile epsilon and omits density/novelty filters too')
 }
 h<-read.csv(file.path(P,'Scientific Reports Manuscript/github_results/summary_results/negative_control_headline.csv'));hh<-list()
 for(method in unique(orig$method))for(m in c('AnyMicroDiscovery','SpuriousExtraClusters','MajorARI')){a<-orig[orig$method==method,];v<-mean(a[[m]],na.rm=TRUE);target<-h[[m]][h$method==method];stopifnot(length(target)==1,abs(v-target)<1e-11);hh[[length(hh)+1]]<-data.frame(method=method,metric=m,n_finite=sum(is.finite(a[[m]])),reproduced=v,submitted=target,absolute_difference=abs(v-target))}
 wr(bind(hh),'original_headline_checks.csv','evidence')
 wr(ocells,'original_negative_design.csv');wr(bind(os),'original_negative_cell_discoveries.csv');wr(bind(opool),'original_negative_pooled.csv');wr(bind(omatch),'original_paired_ablation.csv')
 cp<-binom.test(3,10000)$conf.int
 cl<-vapply(split(o,o$cell_id),function(a)binom.test(sum(a$AnyMicroDiscovery),nrow(a),conf.level=1-.05/10)$conf.int[1],numeric(1));cu<-vapply(split(o,o$cell_id),function(a)binom.test(sum(a$AnyMicroDiscovery),nrow(a),conf.level=1-.05/10)$conf.int[2],numeric(1))
 wr(data.frame(n_records=10000,n_cells=10,events=3,frequency=.0003,conventional_binomial95_lower=cp[1],conventional_binomial95_upper=cp[2],fixed_equal_cell_mean_conservative95_lower=mean(cl),fixed_equal_cell_mean_conservative95_upper=mean(cu),interpretation='Fixed10-cell target:average simultaneous99.5% per-cell exact intervals gives conservative>=95% bounds by union bound;not a universal noise guarantee'),'original_interval_scope.csv')
 # A-END reconciliation uses nearest exact DGP fields EXCEPT sample size;never a paired contrast.
 reconciliation<-list()
 for(cid in unique(x$cell_id[x$experiment=='A-END'&x$rule=='legacy_rule'])){
  a<-x[x$cell_id==cid,];meta<-a[1,c('scenario','noise_type','p','requested_alpha','n','n_residual','cell_id')]
  candidates<-ocells[ocells$noise_type==a$noise_type[1]&ocells$p==a$p[1]&ocells$requested_alpha==a$requested_alpha[1]&ocells$main_distribution=='normal',]
  stopifnot(nrow(candidates)<=1)
  if(nrow(candidates)){
   ob<-o[o$cell_id==candidates$cell_id,];oldrate<-mean(ob$AnyMicroDiscovery);oid<-candidates$cell_id;on<-nrow(ob);oe<-sum(ob$AnyMicroDiscovery)
  } else{oldrate<-oid<-on<-oe<-NA_real_}
  reconciliation[[length(reconciliation)+1]]<-cbind(meta,tier_records=nrow(a),tier_events=sum(a$AnyMicroDiscovery),tier_frequency=mean(a$AnyMicroDiscovery),original_cell=oid,original_records=on,original_events=oe,original_frequency=oldrate,original_n=if(nrow(candidates))candidates$n else NA,original_residual_size_engine=if(nrow(candidates))candidates$n-floor((1-candidates$requested_alpha)*candidates$n) else NA,comparison=if(nrow(candidates))'Same named noise-generator formula,p,trimming,K2,normal elliptical macros;different total n1000 versus500,residual sizes,seeds and independent fits. Descriptive comparison only.' else 'No original cell matching this noise mechanism,dimension and trimming;not treated as a matched population',source='Original negative_control_replications.rds plus Tier1 A-END;separate populations')
 }
 wr(bind(reconciliation),'original_vs_tier1_reconciliation.csv')
 availability<-unique(x[c('experiment','mode','noise_type','rule')]);availability$n_cells<-vapply(seq_len(nrow(availability)),function(i)sum(vapply(split(x,x$cell_id),function(a)all(vapply(names(availability)[1:4],function(nm)a[[nm]][1]==availability[[nm]][i],logical(1))),logical(1))),integer(1));wr(availability,'rule_noise_availability.csv')
 saveRDS(list(negative=x,original_validated=o,original_unvalidated=u),file.path(C,'evidence/audited_scalar_records.rds'),compress='xz')
 wr(data.frame(check=c('235 immutable inputs','Frozen Tier1 design and full results','All9880 task blocks represented','All design fields and deterministic seeds','63800 TSR no-minority records and integer assignment counts','36 verified paired rule comparisons','Original80000 rows and66 export fields','Original3/10000 and10000/10000 discovery events','Original native-event identification from strict max rule','No fitting or simulation'),status='PASS',detail=c('SHA256 checked before audit','282 configurations;98800 OK records','tasks_tier1.csv checked;cell/rep keys complete','Every common field matches;seed formulas reproduce','All QHat==Spurious;native count/event agree;M/noise denominators explicit',paste(counter,'contrasts;all nonrule inputs and macro diagnostics identical'),'Blank serialized text normalized only for comparison','All original successful rows retained;ten fixed cells kept separate','Exactly3 derived native events,identical to external events;accepted sizes remain unavailable','Read-only summaries and deterministic calculations only')),'analysis_checks.csv','evidence')
 writeLines(capture.output(sessionInfo()),file.path(C,'evidence/R_session_info.txt'))
 cat('R2C1_READ_ONLY_AUDIT_OK',nrow(x),'negative TSR rows;',counter,'paired contrasts;NO_FITS\n')
 print(bind(overview)[bind(overview)$experiment%in%c('A-END','A-COND'),c('scenario','rule','external_events','native_events','native_rate','mean_residual_assigned')],row.names=FALSE)
 print(bind(os)[bind(os)$method=='TSRClust_validated',c('cell_id','noise_type','p','requested_alpha','events')],row.names=FALSE)
}
main()
