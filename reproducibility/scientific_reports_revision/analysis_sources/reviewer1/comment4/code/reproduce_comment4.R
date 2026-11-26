# Reviewer1 Comment4 read-only audit: no generated observations, fits or installs.
stopifnot(as.character(getRversion())=='4.5.2')
f<-grep('^--file=',commandArgs(),value=TRUE)[1];code<-dirname(normalizePath(sub('^--file=','',f)));c4<-dirname(code);revision<-dirname(dirname(c4));project<-dirname(revision)
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
inv<-read.csv(file.path(c4,'evidence/source_inventory.csv'));for(i in seq_len(nrow(inv))){p<-if(nzchar(inv$project_relative[i]))file.path(project,inv$project_relative[i]) else archive;stopifnot(sha(p)==inv$sha256[i])}
write<-function(x,n,where='tables')write.csv(x,file.path(c4,where,n),row.names=FALSE,na='NA')
bind<-function(xs)do.call(rbind,xs)
mean_ci<-function(v,binary=FALSE){
 v<-v[is.finite(v)];n<-length(v)
 if(!n)return(c(n_finite=0,events=NA,mean=NA,sd=NA,mcse=NA,lower=NA,upper=NA,median=NA))
 se<-if(n>1)sd(v)/sqrt(n) else NA_real_;ci<-if(n>1)mean(v)+c(-1,1)*qt(.975,n-1)*se else c(NA,NA)
 if(binary){stopifnot(all(v%in%c(0,1)));ci<-binom.test(sum(v),n)$conf.int}
 c(n_finite=n,events=if(binary)sum(v) else NA,mean=mean(v),sd=if(n>1)sd(v) else NA,mcse=se,lower=ci[1],upper=ci[2],median=median(v))
}
stratified_ci<-function(v,cell){
 ok<-is.finite(v);v<-v[ok];cell<-cell[ok];n<-length(v);z<-mean_ci(v)
 if(!n)return(z)
 groups<-split(v,cell);nc<-lengths(groups);vv<-vapply(groups,function(a)if(length(a)>1)length(a)*var(a)/n^2 else NA_real_,numeric(1))
 se<-if(all(is.finite(vv)))sqrt(sum(vv)) else NA_real_
 df<-if(is.finite(se)&&se>0)sum(vv)^2/sum(vv^2/(nc-1)) else NA_real_
 ci<-if(is.finite(se)&&se==0)rep(mean(v),2) else if(is.finite(df))mean(v)+c(-1,1)*qt(.975,df)*se else c(NA,NA)
 z[c('mcse','lower','upper')]<-c(se,ci);z
}
paired_ci<-function(a,b,binary=FALSE){
 ok<-is.finite(a)&is.finite(b);a<-a[ok];b<-b[ok];n<-length(a);d<-a-b
 z<-mean_ci(d);z<-z[c('n_finite','mean','sd','mcse','lower','upper')]
 up<-down<-NA_integer_
 if(binary&&n){
  stopifnot(all(a%in%c(0,1)),all(b%in%c(0,1)));up<-sum(a==1&b==0);down<-sum(a==0&b==1)
  # Two 97.5% exact binomial margins give simultaneous >=95% discordance bounds.
  cu<-binom.test(up,n,conf.level=.975)$conf.int;cd<-binom.test(down,n,conf.level=.975)$conf.int
  z[c('lower','upper')]<-c(cu[1]-cd[2],cu[2]-cd[1])
 }
 c(z,n_up=up,n_down=down)
}
base_metrics<-c('MajorARI','GlobalARI','MajorRetention','MinorityClusterF1','MinorityDetectionF1','MinorityDetectionF1_all_extra_v2','MacroFalsePositiveReassignment','FalseExtraAssignmentRate','AnyMicroDiscovery','SpuriousExtraClusters','QHat','any_native_micro','n_accepted_micro','end_to_end_seconds','macro_seconds','fit_seconds')
binary_metrics<-c('AnyMicroDiscovery','any_native_micro','AnySpuriousExtra','K_correct_macro','K_above_macro','K_below_macro','K_correct_total','K_above_total','K_at_upper_boundary')
work<-tempfile('c4_audit_');dir.create(work)
run_audit<-function(){
 on.exit(unlink(work,recursive=TRUE));untar(archive,files='replications.csv',exdir=work);all<-read.csv(file.path(work,'replications.csv'))
 stopifnot(nrow(all)==98800,all(all$status=='OK'),length(unique(all$cell_id))==282)
 design<-read.csv(file.path(project,'revision_scientific_reports/config/design.csv'));manifest<-design[design$tier==1,]
 stopifnot(nrow(manifest)==282,!anyDuplicated(all[c('cell_id','replication')]))
 # Validate every stored design field against the frozen design, not family names.
 fields<-intersect(names(design),names(all));mm<-match(all$cell_id,design$cell_id)
 for(nm in fields){a<-all[[nm]];b<-design[[nm]][mm];stopifnot(identical(is.na(a),is.na(b)),all(as.character(a)==as.character(b),na.rm=TRUE))}
 e<-new.env();sys.source(file.path(project,'revision_scientific_reports/R/common.R'),e)
 stopifnot(all(all$data_seed==mapply(e$seed_id,all$dataset_id,all$replication,MoreArgs=list(stream=0L))),all(all$algorithm_seed==mapply(e$seed_id,all$dataset_id,all$replication,all$algorithm_repeat,MoreArgs=list(stream=1L))),all(all$null_seed==mapply(e$seed_id,all$dataset_id,all$replication,all$null_repeat,MoreArgs=list(stream=2L))))
 family<-list()
 for(exp in unique(all$experiment))for(m in unique(all$method[all$experiment==exp]))for(info in unique(all$K_information[all$experiment==exp&all$method==m])){
  a<-all[all$experiment==exp&all$method==m&all$K_information==info,]
  family[[length(family)+1]]<-data.frame(experiment=exp,method=m,K_information=info,regime=paste(unique(a$regime),collapse=';'),nominal_configurations=length(unique(a$cell_id)),n_records=nrow(a),unique_dataset_repetitions=nrow(unique(a[c('dataset_id','replication')])),kmax_values=paste(sort(unique(a$kmax)),collapse=';'),p_values=paste(sort(unique(a$p)),collapse=';'),noise_types=paste(sort(unique(a$noise_type)),collapse=';'),scenario_values=paste(unique(a$scenario),collapse=';'),source='Tier1 archive replications.csv')
 }
 write(bind(family),'tier1_family_inventory.csv')
 x<-all[all$experiment%in%c('D1-E','D2-E'),];stopifnot(nrow(x)==5600)
 x$AnySpuriousExtra<-as.numeric(x$SpuriousExtraClusters>0)
 # Meaningful only for selectors. Density methods' estimated_K is the unchanged
 # initial variable in the worker, NOT a measured number of their clusters.
 sel<-x$experiment=='D2-E';stopifnot(all(x$estimated_K[sel]%in%1:7),all(x$estimated_K[sel]<=x$kmax[sel]))
 for(nm in setdiff(binary_metrics,c('AnyMicroDiscovery','any_native_micro','AnySpuriousExtra')))x[[nm]]<-NA_real_
 x$K_correct_macro[sel]<-as.numeric(x$estimated_K[sel]==x$main_k[sel]);x$K_above_macro[sel]<-as.numeric(x$estimated_K[sel]>x$main_k[sel]);x$K_below_macro[sel]<-as.numeric(x$estimated_K[sel]<x$main_k[sel]);x$K_correct_total[sel]<-as.numeric(x$estimated_K[sel]==x$main_k[sel]+x$q[sel]);x$K_above_total[sel]<-as.numeric(x$estimated_K[sel]>x$main_k[sel]+x$q[sel]);x$K_at_upper_boundary[sel]<-as.numeric(x$estimated_K[sel]==x$kmax[sel])
 metrics<-c(base_metrics,binary_metrics[!binary_metrics%in%base_metrics])
 summaries<-list();distribution<-list();conditional<-list();rates<-list()
 for(cid in unique(x$cell_id)){
  a<-x[x$cell_id==cid,];n<-nrow(a);stopifnot(n==100)
  meta<-data.frame(experiment=a$experiment[1],scenario=a$scenario[1],method=a$method[1],K_information=a$K_information[1],regime=a$regime[1],cell_id=cid,dataset_id=a$dataset_id[1],main_k=a$main_k[1],q=a$q[1],kmax=a$kmax[1],requested_alpha=a$requested_alpha[1],n_independent_configurations=1,n_records=n,n_unique_datasets=n,paired_inference=FALSE,source='Tier1 archive replications.csv')
  for(m in metrics)summaries[[length(summaries)+1]]<-cbind(meta,metric=m,as.data.frame(as.list(mean_ci(a[[m]],m%in%binary_metrics))),uncertainty=if(m%in%binary_metrics)'Exact binomial 95% descriptive cell interval' else 'Pointwise Monte Carlo t interval; one fixed cell')
  if(a$experiment[1]=='D2-E'){
   for(k in seq_len(a$kmax[1]))distribution[[length(distribution)+1]]<-cbind(meta,selected_K=k,count=sum(a$estimated_K==k),frequency=mean(a$estimated_K==k),frequency_lower=binom.test(sum(a$estimated_K==k),n)$conf.int[1],frequency_upper=binom.test(sum(a$estimated_K==k),n)$conf.int[2],uncertainty='Pointwise exact binomial95%; category intervals not simultaneous',meaning=if(a$method[1]=='TSRClust')'Stage-I selected components; target dominant K' else 'Full-data t components; total groups K+Q is only an external reference')
   for(m in c('estimated_K',setdiff(binary_metrics,c('AnyMicroDiscovery','any_native_micro','AnySpuriousExtra'))))rates[[length(rates)+1]]<-cbind(meta,metric=m,as.data.frame(as.list(mean_ci(a[[m]],m%in%binary_metrics))))
   for(k in sort(unique(a$estimated_K)))for(m in c('MajorARI','MinorityClusterF1','AnyMicroDiscovery','any_native_micro','SpuriousExtraClusters'))conditional[[length(conditional)+1]]<-cbind(meta,selected_K=k,metric=m,n_conditional=sum(a$estimated_K==k),as.data.frame(as.list(mean_ci(a[[m]][a$estimated_K==k],m%in%binary_metrics))),caution='Post-selection descriptive subset, not causal or independent validation')
  }
 }
 write(bind(summaries),'performance_summary.csv');write(bind(distribution),'selected_K_distribution.csv');write(bind(rates),'K_selection_rates.csv');write(bind(conditional),'conditional_selection_performance.csv')
 # Every contrast receives an explicit pairing proof over ALL DGP fields and RNG IDs.
 dgp<-c('n','p','main_k','separation','main_structure','main_distribution','main_weights','residual_total','minority_ratio','q','minority_shape','minority_sd','minority_proximity','noise_type','mode','nR')
 streams<-c('data_seed','algorithm_seed','null_seed','algorithm_repeat','null_repeat','algorithm_version','metric_version')
 pairout<-list();proof<-list();counter<-0L
 compare<-function(a,b,label,family,selected_metrics=metrics){
  stopifnot(nrow(a)==100,nrow(b)==100);a<-a[order(a$replication),];b<-b[order(b$replication),]
  stopifnot(identical(a$replication,b$replication),all(a$dataset_id==b$dataset_id),all(vapply(c(dgp,streams),function(z)identical(a[[z]],b[[z]]),logical(1))))
  counter<<-counter+1L;id<-sprintf('C4-%03d',counter)
  proof[[length(proof)+1]]<<-data.frame(contrast_id=id,family=family,label=label,scenario=a$scenario[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],new_method=a$method[1],baseline_method=b$method[1],new_kmax=a$kmax[1],baseline_kmax=b$kmax[1],n_pairs=100,unique_dataset_ids=1,dataset_id=a$dataset_id[1],all_DGP_fields_equal=TRUE,data_seed_equal=TRUE,algorithm_start_seed_equal=TRUE,null_start_seed_equal=TRUE,repeat_IDs_equal=TRUE,scientific_version_equal=TRUE,pairing='Exact under frozen generator/DGP/seed contract; actual X hashes not exported',random_draw_caveat='Same initial seeds, not identical realized streams after model-search branches',source='Tier1 archive replications.csv; frozen design and generator source')
  for(m in selected_metrics){z<-paired_ci(a[[m]],b[[m]],m%in%binary_metrics)
   pairout[[length(pairout)+1]]<<-cbind(data.frame(contrast_id=id,family=family,label=label,scenario=a$scenario[1],new_experiment=a$experiment[1],baseline_experiment=b$experiment[1],new_method=a$method[1],baseline_method=b$method[1],new_information=a$K_information[1],baseline_information=b$K_information[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],new_kmax=a$kmax[1],baseline_kmax=b$kmax[1],n_independent_configurations_per_arm=1,n_records_per_arm=100,n_pairs=100,metric=m,mean_new=if(any(is.finite(a[[m]])))mean(a[[m]],na.rm=TRUE) else NA,mean_baseline=if(any(is.finite(b[[m]])))mean(b[[m]],na.rm=TRUE) else NA),as.data.frame(as.list(z)),direction='New minus baseline',paired=TRUE,uncertainty=if(m%in%binary_metrics)'Conservative exact paired-discordance 95% bound; two97.5% binomial margins' else 'Pointwise paired t Monte Carlo95%; degenerate variance not a population guarantee',source='Tier1 archive replications.csv')
  }
 }
 for(sc in unique(x$scenario)){
  d1<-x[x$experiment=='D1-E'&x$scenario==sc,];d2<-x[x$experiment=='D2-E'&x$scenario==sc,]
  known<-d1[d1$method=='TSRClust',];fixedt<-d1[d1$method=='tEIGEN_fixedK',]
  for(comp in setdiff(unique(d1$method),'TSRClust'))compare(known,d1[d1$method==comp,],paste('TSRClust minus',comp),'Setting A controlled inputs',base_metrics)
  for(km in c(3,5,7)){
   a<-d2[d2$method=='TSRClust'&d2$kmax==km,];b<-d2[d2$method=='tEIGEN_BIC'&d2$kmax==km,]
   compare(a,b,'Selected TSRClust minus BIC tEIGEN','Setting B common K/G range',base_metrics)
   compare(a,known,'Selected TSRClust minus true-K TSRClust','Setting C full selection-pipeline effect',base_metrics)
   compare(b,fixedt,'BIC tEIGEN minus fixed-G tEIGEN','Setting C tEIGEN component-search effect',base_metrics)
  }
  for(m in unique(d2$method))for(km in c(5,7))compare(d2[d2$method==m&d2$kmax==km,],d2[d2$method==m&d2$kmax==3,],'Wider range minus range1:3','Setting C range sensitivity',c(metrics,'estimated_K'))
 }
 write(bind(pairout),'paired_performance_differences.csv');write(bind(proof),'pairing_checks.csv','evidence')
 knownsel<-bind(pairout);write(knownsel[knownsel$family=='Setting C full selection-pipeline effect',],'known_vs_selected_K.csv')
 write(knownsel[knownsel$family=='Setting C range sensitivity',],'range_sensitivity.csv')
 # Independence audit: range sweeps re-use these data; do NOT multiply n by grids/methods.
 seen<-unique(x[c('dataset_id','replication','data_seed')]);stopifnot(nrow(seen)==400)
 write(data.frame(experiments='D1-E+D2-E',n_scalar_records=nrow(x),n_nominal_cells=length(unique(x$cell_id)),n_unique_generated_datasets=nrow(seen),n_scenarios=length(unique(x$scenario)),replications_per_scenario=100,warning='56 design cells share400 datasets; no independent pooling across methods or search ranges'),'independence_audit.csv','evidence')
 # Original submitted method information checks and selected-G distributions.
 origchecks<-list();origdist<-list();origperf<-list();crossdesign<-list()
 for(st in c('architecture','stress','negative_control')){
  a<-readRDS(file.path(project,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined',paste0(st,'_replications.rds')))
  # Original/Tier1 cells can share DGP settings yet have disjoint seed namespaces.
  for(sc in unique(x$scenario)){
   target<-x[x$scenario==sc,][1,];candidates<-a[a$method=='TSRClust_validated',]
   for(cid in unique(candidates$cell_id)){
    zc<-candidates[candidates$cell_id==cid,]
    equal<-all(vapply(setdiff(dgp,c('mode','nR')),function(nm)as.character(zc[[nm]][1])==as.character(target[[nm]]),logical(1)))
    if(equal){seed_overlap<-length(intersect(zc$seed,unique(x$data_seed[x$scenario==sc])));stopifnot(seed_overlap==0)
     crossdesign[[length(crossdesign)+1]]<-data.frame(original_study=st,original_cell=cid,tier1_scenario=sc,tier1_dataset_id=target$dataset_id,n_original_datasets=nrow(zc),n_tier1_datasets=100,n_shared_data_seeds=seed_overlap,comparison='Matched DGP design but distinct data seeds; not paired; descriptive only')
    }
   }
  }
  for(m in unique(a$method)){
   z<-a[a$method==m,];good<-z$status=='OK'
   if(m%in%c('TCLUST','TSRClust_validated','TSRClust_unvalidated','Trimmed_kmeans','OTRIMLE_fixedK'))stopifnot(all(z$requested_k==z$main_k),all(z$requested_alpha==z$residual_total))
   if(m=='Oracle_TCLUST')stopifnot(all(z$fitted_k[good]==z$main_k[good]+z$q[good]),all(abs(z$fitted_alpha[good]-(z$residual_total[good]-z$minority_ratio[good]))<1e-12))
   origchecks[[length(origchecks)+1]]<-data.frame(study=st,method=m,n_records=nrow(z),n_OK=sum(good),n_configurations=length(unique(z$cell_id)),requested_k_equals_true_macro=sum(z$requested_k==z$main_k,na.rm=TRUE),requested_alpha_equals_true_residual=sum(z$requested_alpha==z$residual_total,na.rm=TRUE),n_fitted_K_finite=sum(is.finite(z$fitted_k)),n_selected_G_finite=sum(is.finite(z$selected_g)),n_E2E_runtime_finite=sum(is.finite(z$end_to_end_seconds)),interpretation='Requested fields are design metadata for every method; not necessarily actual consumed arguments. Read method information audit.')
   for(co in unique(if(st=='architecture')ifelse(z$q>0,'M1-M4','N0') else st)){
    b<-if(st=='architecture')z[if(co=='N0')z$q==0 else z$q>0,] else z
    for(metric in intersect(base_metrics,names(b)))origperf[[length(origperf)+1]]<-cbind(study=st,cohort=co,method=m,n_records=nrow(b),metric=metric,as.data.frame(as.list(stratified_ci(b[[metric]],b$cell_id))),uncertainty='Fixed-cell stratified Monte Carlo95%; Welch df; original cohort kept separate')
   }
   if(m=='tEIGEN_BIC'){
    gmax<-pmin(6,z$main_k+pmax(z$q,1)+2);stopifnot(all(z$selected_g[good]>=1),all(z$selected_g[good]<=gmax[good]))
    for(co in unique(paste(z$main_k,z$q,gmax,sep='/'))){b<-z[paste(z$main_k,z$q,gmax,sep='/')==co,];for(g in sort(unique(b$selected_g[is.finite(b$selected_g)])))origdist[[length(origdist)+1]]<-data.frame(study=st,method=m,main_k=b$main_k[1],q=b$q[1],gmax=min(6,b$main_k[1]+max(b$q[1],1)+2),selected_G=g,n_records=nrow(b),n_finite=sum(is.finite(b$selected_g)),count=sum(b$selected_g==g,na.rm=TRUE),frequency=mean(b$selected_g==g,na.rm=TRUE))}
   }
  }
 }
 if(length(crossdesign))write(bind(crossdesign),'cross_study_design_matches.csv')
 write(bind(origchecks),'original_information_checks.csv');write(bind(origdist),'original_tEIGEN_G_distribution.csv');write(bind(origperf),'original_performance_descriptive.csv')
 # Existing old automatic run inspected separately; never evidence for current D2.
 old<-readRDS(file.path(project,'TSRClust_SciReports_TRUBA/results/combined/automacro_replications.rds'))
 write(data.frame(path='TSRClust_SciReports_TRUBA/results/combined/automacro_replications.rds',n_rows=nrow(old),n_fields=ncol(old),methods=paste(unique(old$method),collapse=';'),p=paste(unique(old$p),collapse=';'),n=paste(unique(old$n),collapse=';'),has_end_to_end_runtime='end_to_end_seconds'%in%names(old),status='Historical auxiliary record; excluded from current benchmark performance/paired inference'),'historical_auto_availability.csv','evidence')
 write(data.frame(check=c(paste('all',nrow(inv),'immutable scientific input hashes'),'Tier1 dimensions and frozen design fields','all data/algorithm/null seed formulas','all contrast DGP/RNG pairing checks','selected K within declared ranges','400 independent datasets across56 D/E cells','original fixed-K/alpha and Oracle inputs','original tEIGEN truth-assisted candidate-range formula'),status='PASS',detail=c('SHA256 verified before any calculation','98800 OK records;282 cells;5600 D/E rows','Source seed_id matches every record',paste(counter,'explicit contrast proofs'),'All2400 D2 selected values finite/in range','Range/method repeats not pooled as independent','Actual fitted/used inputs audited separately from row metadata','All original finite selected G within computed range')),'validation_checks.csv','evidence')
 # Regression fixtures exercise uncertainty logic, not model fitting.
 z<-paired_ci(c(0,1,1,0),c(0,0,1,1),TRUE);stopifnot(z['mean']==0,z['n_up']==1,z['n_down']==1,z['lower']<0,z['upper']>0)
 z<-paired_ci(rep(0,100),rep(0,100),TRUE);stopifnot(z['mean']==0,z['lower']<0,z['upper']>0)
 v<-matrix(c(5,1,1,3),nrow=2);pos<-which(v==min(v),arr.ind=TRUE)[1,];stopifnot(identical(as.integer(pos),c(2L,1L)),which.min(c(1,1))==1L)
 write(data.frame(test=c('exact paired discordance nondegenerate zero events','first column-major IC tie, then first alpha candidate'),status='PASS'),'fixed_rule_regressions.csv','evidence')
 writeLines(capture.output(sessionInfo()),file.path(c4,'evidence/R_session_info.txt'))
 cat('COMMENT4_READONLY_AUDIT_OK',counter,'paired contrasts;',nrow(bind(pairout)),'metric contrast rows\n')
}
run_audit()
