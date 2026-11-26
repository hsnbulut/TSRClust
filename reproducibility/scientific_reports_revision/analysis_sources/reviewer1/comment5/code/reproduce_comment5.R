# Read-only Reviewer 1 Comment 5 audit. No fitting, simulation or installation.
stopifnot(as.character(getRversion())=='4.5.2')
f<-grep('^--file=',commandArgs(),value=TRUE)[1]
code<-dirname(normalizePath(sub('^--file=','',f)));C<-dirname(code);R<-dirname(dirname(C));P<-dirname(R)
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
inv<-read.csv(file.path(C,'evidence/source_inventory.csv'))
for(i in seq_len(nrow(inv))){p<-if(nzchar(inv$project_relative[i]))file.path(P,inv$project_relative[i]) else archive;stopifnot(sha(p)==inv$sha256[i])}
wr<-function(x,n,where='tables')write.csv(x,file.path(C,where,n),row.names=FALSE,na='NA')
bind<-function(xs)do.call(rbind,xs)
ci<-function(v,binary=FALSE){
 v<-v[is.finite(v)];n<-length(v)
 if(!n)return(c(n_finite=0,events=NA,mean=NA,sd=NA,mcse=NA,lower=NA,upper=NA,median=NA))
 se<-if(n>1)sd(v)/sqrt(n) else NA_real_;b<-if(n>1)mean(v)+c(-1,1)*qt(.975,n-1)*se else c(NA,NA)
 if(binary){stopifnot(all(v%in%c(0,1)));b<-binom.test(sum(v),n)$conf.int}
 c(n_finite=n,events=if(binary)sum(v) else NA,mean=mean(v),sd=if(n>1)sd(v) else NA,mcse=se,lower=b[1],upper=b[2],median=median(v))
}
paired<-function(a,b,binary=FALSE){
 ok<-is.finite(a)&is.finite(b);a<-a[ok];b<-b[ok];z<-ci(a-b)[c('n_finite','mean','sd','mcse','lower','upper')];up<-down<-NA_real_
 if(binary&&length(a)){stopifnot(all(a%in%c(0,1)),all(b%in%c(0,1)));up<-sum(a==1&b==0);down<-sum(a==0&b==1);cu<-binom.test(up,length(a),conf.level=.975)$conf.int;cd<-binom.test(down,length(a),conf.level=.975)$conf.int;z[c('lower','upper')]<-c(cu[1]-cd[2],cu[2]-cd[1])}
 c(z,n_up=up,n_down=down)
}
metrics<-c('MinorityClusterF1','MinorityDetectionF1_all_extra_v2','AnyMicroDiscovery','any_native_micro','SpuriousExtraClusters','FalseExtraAssignmentRate','MacroFalsePositiveReassignment','NoiseRecall','MajorARI','n_accepted_micro','micro_seconds','end_to_end_seconds')
binary<-c('AnyMicroDiscovery','any_native_micro')
work<-tempfile('c5_readonly_');dir.create(work)
audit<-function(){
 on.exit(unlink(work,recursive=TRUE));untar(archive,files='replications.csv',exdir=work);all<-read.csv(file.path(work,'replications.csv'))
 design<-read.csv(file.path(P,'revision_scientific_reports/config/design.csv'));stopifnot(nrow(all)==98800,all(all$status=='OK'),length(unique(all$cell_id))==282,!anyDuplicated(all[c('cell_id','replication')]))
 mm<-match(all$cell_id,design$cell_id)
 for(nm in intersect(names(design),names(all))){a<-all[[nm]];b<-design[[nm]][mm];stopifnot(identical(is.na(a),is.na(b)),all(as.character(a)==as.character(b),na.rm=TRUE))}
 e<-new.env();sys.source(file.path(P,'revision_scientific_reports/R/common.R'),e)
 stopifnot(all(all$data_seed==mapply(e$seed_id,all$dataset_id,all$replication,MoreArgs=list(stream=0L))),all(all$algorithm_seed==mapply(e$seed_id,all$dataset_id,all$replication,all$algorithm_repeat,MoreArgs=list(stream=1L))),all(all$null_seed==mapply(e$seed_id,all$dataset_id,all$replication,all$null_repeat,MoreArgs=list(stream=2L))))
 x<-all[all$method=='TSRClust',];x$minPts<-pmax(5,x$minpts_multiplier*x$p)
 stopifnot(nrow(x)==94800,all(is.finite(x$n_residual)),all(x$n_residual>=x$minPts),all(x$rejection_reason!='infeasible_minPts'))
 x$implied_trimmed<-ifelse(x$mode=='residual_only',NA,x$n-floor((1-x$requested_alpha)*x$n))
 stopifnot(all(x$n_residual==x$implied_trimmed,na.rm=TRUE))
 inventory<-list();reasons<-list()
 for(cid in unique(x$cell_id)){
  a<-x[x$cell_id==cid,];meta<-data.frame(experiment=a$experiment[1],scenario=a$scenario[1],cell_id=cid,dataset_id=a$dataset_id[1],mode=a$mode[1],p=a$p[1],n=a$n[1],q=a$q[1],minority_ratio=a$minority_ratio[1],population_minority_total=if(a$mode[1]=='residual_only')NA else round(a$n[1]*a$minority_ratio[1]),requested_alpha=a$requested_alpha[1],minpts_multiplier=a$minpts_multiplier[1],minPts=a$minPts[1],density=a$density[1],n_records=nrow(a),n_residual_finite=sum(is.finite(a$n_residual)),n_residual_min=min(a$n_residual),n_residual_max=max(a$n_residual),n_below=sum(a$n_residual<a$minPts),n_equal=sum(a$n_residual==a$minPts),n_above=sum(a$n_residual>a$minPts),n_native_accepted=sum(a$any_native_micro),source='Tier1 replications.csv; actual n_residual; population U from frozen DGP, not residual truth counts')
  inventory[[length(inventory)+1]]<-meta
  for(why in unique(a$rejection_reason))reasons[[length(reasons)+1]]<-cbind(meta,rejection_reason=why,count=sum(a$rejection_reason==why),fraction=mean(a$rejection_reason==why))
 }
 wr(bind(inventory),'tier1_residual_inventory.csv');wr(bind(reasons),'tier1_gate_outcomes.csv')
 wr(data.frame(field=names(all),available=TRUE,source='Tier1 replications.csv header'),'export_schema.csv','evidence')
 fields<-list(
  c('p,n,requested_alpha,minpts_multiplier','DIRECT','DGP/fit configuration; actual minPts is max(5,multiplier*p)'),
  c('n_residual','DIRECT','Realized Stage-I zero-label count (residual-only: all supplied residual points)'),
  c('population_minority_total','DERIVED','round(n*minority_ratio) from DGP; no residual truth composition'),
  c('early_stop','DERIVED','Compare stored n_residual with actual minPts; no infeasible_minPts observed'),
  c('DBSCAN/core/candidate counts','UNAVAILABLE','scan is returned in memory but not serialized in scalar archive'),
  c('size/density/novelty pass counts','UNAVAILABLE','no_geometric_candidate aggregates all these and empty-grid/no-cluster cases'),
  c('global_p,rejection_reason','DIRECT','Final pipeline reason and scalar global p; not individual candidate evidence'),
  c('null simulation executed','PARTLY IDENTIFIABLE','no_geometric_candidate/infeasible_minPts return before null loop; other reasons after null loop'),
  c('any_native_micro,n_accepted_micro','DIRECT','Native Stage-II event/count; distinct from external AnyMicroDiscovery'),
  c('MinorityClusterF1,MinorityDetectionF1_all_extra_v2,MajorARI','DIRECT','Metric-specific finite masks; negative-control minority F1 undefined'),
  c('micro_seconds,end_to_end_seconds','DIRECT','Measured fitting time; not SLURM allocation or a hardware-independent speed claim'),
  c('minority/core memberships,realized epsilon grid','UNAVAILABLE','Do not infer from scalar scores, population mass or result-object byte count'))
 wr(bind(lapply(fields,function(v)data.frame(quantity=v[1],availability=v[2],interpretation=v[3]))),'diagnostic_availability.csv')
 # C main effect: only minPts multiplier varies. Interactions are kept separate.
 z<-x[x$experiment=='C' & x$reference=='legacy_empirical' & x$rule=='legacy_rule' & x$novelty==.995 & x$B==49 & x$density%in%c(1.5,2,3),]
 z<-z[(z$density==2 & z$minpts_multiplier%in%c(1,2,3,4)) | (z$density!=2 & z$minpts_multiplier%in%c(1,4)),]
 stopifnot(nrow(z)==6400,length(unique(z$cell_id))==32)
 sums<-list();pout<-list();proof<-list();counter<-0L
 for(cid in unique(z$cell_id)){
  a<-z[z$cell_id==cid,];stopifnot(nrow(a)==200)
  meta<-data.frame(family=if(a$density[1]==2)'Controlled default-density multiplier sweep' else 'Separate fixed-density interaction slice',scenario=a$scenario[1],cell_id=cid,dataset_id=a$dataset_id[1],p=a$p[1],n=a$n[1],requested_alpha=a$requested_alpha[1],n_residual_min=min(a$n_residual),n_residual_max=max(a$n_residual),minpts_multiplier=a$minpts_multiplier[1],minPts=a$minPts[1],density=a$density[1],novelty=.995,B=49,n_records=nrow(a),n_unique_generated_datasets=nrow(a),source='Tier1 experiment C; fixed macro count; uniform generated noise')
  for(m in metrics)sums[[length(sums)+1]]<-cbind(meta,metric=m,as.data.frame(as.list(ci(a[[m]],m%in%binary))),uncertainty=if(m%in%binary)'Pointwise exact binomial95%' else 'Pointwise Monte Carlo t95%; metric-specific finite mask')
 }
 equal_fields<-setdiff(intersect(names(design),names(z)),c('cell_id','minpts_multiplier'))
 equal_fields<-union(equal_fields,c('dataset_id','replication','data_seed','algorithm_seed','null_seed','algorithm_version','metric_version','n_residual','macro_objective','macro_code'))
 compare<-function(a,b,family){
  a<-a[order(a$replication),];b<-b[order(b$replication),];stopifnot(nrow(a)==200,nrow(b)==200)
  for(nm in equal_fields)stopifnot(identical(a[[nm]],b[[nm]]))
  counter<<-counter+1L;id<-sprintf('C5-%02d',counter)
  proof[[length(proof)+1]]<<-data.frame(contrast_id=id,family=family,scenario=a$scenario[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],new_multiplier=a$minpts_multiplier[1],baseline_multiplier=b$minpts_multiplier[1],density=a$density[1],n_pairs=200,dataset_id=a$dataset_id[1],all_other_design_fields_equal=TRUE,all_three_starting_seeds_equal=TRUE,macro_objective_code_and_R0_equal=TRUE,pairing='Frozen DGP/seed identity; X hashes and realized post-branch draws not exported',source='Tier1 C replications.csv and frozen design')
  for(m in metrics)pout[[length(pout)+1]]<<-cbind(data.frame(contrast_id=id,family=family,scenario=a$scenario[1],new_cell=a$cell_id[1],baseline_cell=b$cell_id[1],new_multiplier=a$minpts_multiplier[1],baseline_multiplier=b$minpts_multiplier[1],new_minPts=a$minPts[1],baseline_minPts=b$minPts[1],density=a$density[1],n_records_per_arm=200,n_pairs=200,metric=m,mean_new=if(any(is.finite(a[[m]])))mean(a[[m]][is.finite(a[[m]])]) else NA,mean_baseline=if(any(is.finite(b[[m]])))mean(b[[m]][is.finite(b[[m]])]) else NA),as.data.frame(as.list(paired(a[[m]],b[[m]],m%in%binary))),direction='New minus baseline',uncertainty=if(m%in%binary)'Conservative exact paired-discordance95%; two97.5% binomial margins' else 'Pointwise paired Monte Carlo t95%; not multiplicity-adjusted',source='Tier1 C replications.csv')
 }
 for(sc in unique(z$scenario)){
  for(mp in c(1,3,4))compare(z[z$scenario==sc&z$density==2&z$minpts_multiplier==mp,],z[z$scenario==sc&z$density==2&z$minpts_multiplier==2,],'Controlled default-density multiplier sweep')
  for(d in c(1.5,3))compare(z[z$scenario==sc&z$density==d&z$minpts_multiplier==4,],z[z$scenario==sc&z$density==d&z$minpts_multiplier==1,],'Separate fixed-density interaction slice')
 }
 stopifnot(counter==20);wr(bind(sums),'minpts_cell_performance.csv');wr(bind(pout),'minpts_paired_differences.csv');wr(bind(proof),'pairing_checks.csv','evidence')
 wr(data.frame(n_main_cells=16,n_main_records=3200,n_main_contrasts=12,n_interaction_cells=16,n_interaction_records=3200,n_interaction_contrasts=8,n_unique_generated_datasets=nrow(unique(z[c('dataset_id','replication','data_seed')])),warning='Cells/alternatives reuse200 datasets per scenario; no pooling as6400 independent datasets'),'independence_audit.csv','evidence')
 # Existing high-dimensional original results: exact original exports, not an older manuscript.
 hd<-list();orig<-list();origmetrics<-list()
 for(st in c('architecture','stress','negative_control')){
  a<-readRDS(file.path(P,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined',paste0(st,'_replications.rds')));a<-a[a$method=='TSRClust_validated',]
  for(p in sort(unique(a$p)))orig[[length(orig)+1]]<-data.frame(study=st,p=p,n_records=sum(a$p==p),source=paste0('Submitted final_sim_combined/',st,'_replications.rds'))
  for(cid in unique(a$cell_id[a$p>=25])){
   b<-a[a$cell_id==cid,];stopifnot(all(b$status=='OK'),all(b$micro_score==0),all(b$micro_p_value==1));m<-max(5,2*b$p[1]);nr<-b$n[1]-floor((1-b$requested_alpha[1])*b$n[1])
   meta<-data.frame(population=paste('Original submitted',st),cell_id=cid,p=b$p[1],n=b$n[1],requested_alpha=b$requested_alpha[1],minpts_multiplier=2,minPts=m,n_records=nrow(b),n_residual_observed=NA,n_residual_implied=nr,residual_count_basis='Engine-implied exact TCLUST2.2-0 hard trimming; scalar RDS has no residual count',n_below=if(nr<m)nrow(b) else 0,n_equal=if(nr==m)nrow(b) else 0,n_above=if(nr>m)nrow(b) else 0,n_no_geometric_candidate=NA,n_score_zero=sum(b$micro_score==0),n_global_p_one=sum(b$micro_p_value==1),n_native_accepted=NA,n_external_discoveries=sum(b$AnyMicroDiscovery),q=b$q[1],population_minority_total=round(b$n[1]*b$minority_ratio[1]),mechanism=if(nr==m)'Outer StageII entered (finite score0,p1); equality size+density identity prevents positive score' else 'Outer StageII entered (finite score0,p1); no positive observed score; individual failed gates not saved',source=paste0('Submitted final_sim_combined/',st,'_replications.rds'))
   hd[[length(hd)+1]]<-meta
   for(metric in intersect(c('MinorityClusterF1','MinorityDetectionF1','AnyMicroDiscovery','MajorARI','micro_score','micro_p_value','end_to_end_seconds'),names(b)))origmetrics[[length(origmetrics)+1]]<-cbind(data.frame(study=st,cell_id=cid,p=b$p[1],n_records=nrow(b),metric=metric,source=meta$source),as.data.frame(as.list(ci(b[[metric]],metric=='AnyMicroDiscovery'))),interpretation='Original submitted high-dimensional scalar summaries; separate from revision Tier1; minority F1 undefined in negative control')
  }
 }
 for(cid in unique(x$cell_id[x$p>=25])){
  b<-x[x$cell_id==cid,];hd[[length(hd)+1]]<-data.frame(population=paste('Tier1',b$experiment[1],b$scenario[1]),cell_id=cid,p=b$p[1],n=b$n[1],requested_alpha=b$requested_alpha[1],minpts_multiplier=b$minpts_multiplier[1],minPts=b$minPts[1],n_records=nrow(b),n_residual_observed=min(b$n_residual),n_residual_implied=b$implied_trimmed[1],residual_count_basis='Observed n_residual in Tier1 scalar export; min/max checked',n_below=sum(b$n_residual<b$minPts),n_equal=sum(b$n_residual==b$minPts),n_above=sum(b$n_residual>b$minPts),n_no_geometric_candidate=sum(b$rejection_reason=='no_geometric_candidate'),n_score_zero=NA,n_global_p_one=sum(b$global_p==1,na.rm=TRUE),n_native_accepted=sum(b$any_native_micro),n_external_discoveries=sum(b$AnyMicroDiscovery),q=b$q[1],population_minority_total=round(b$n[1]*b$minority_ratio[1]),mechanism=if(all(b$n_residual==b$minPts)&&b$density[1]>1)'Equality size+density identity prevents positive score; not outer early stop' else 'No positive geometry candidate; size/density/novelty/core/grid failure counts unavailable',source='Tier1 replications.csv')
 }
 wr(bind(hd),'high_dimensional_failure_diagnostics.csv');wr(bind(orig),'original_dimension_inventory.csv');wr(bind(origmetrics),'original_highdim_metrics.csv')
 # Analytic/count-only boundary fixtures: no coordinates, nearest-neighbour calls or fitting.
 boundary<-list()
 for(p in c(2,5,10,25,50)){
  m<-max(5,2*p)
  for(nr in c(m-1,m,m+1))boundary[[length(boundary)+1]]<-data.frame(p=p,minPts=m,n_residual=nr,relation=if(nr<m)'below' else if(nr==m)'equal' else 'above',stageII_skipped=nr<m,max_possible_self_inclusive_neighbours=nr,DBSCAN_core_possible=if(nr<m)'NO for every epsilon' else 'YES for suitable geometry/epsilon; not guaranteed on actual grid',knn_order=if(nr<m)NA else m-1,validated_acceptance_at_density2=if(nr<=m)'IMPOSSIBLE' else 'NOT GUARANTEED',reason=if(nr<m)'Neighbour-count feasibility: maximum count below minPts; StageI labels unchanged, residuals remain0' else if(nr==m)'Any size-eligible candidate is allR0; background=within distance, density<=1; default2 fails' else 'Core and candidate size may be feasible; grid,density,novelty and MC can still reject',source='Analytic consequences of audited submitted/pinned implementations; not empirical fitted results',n_empirical_fits=0)
 }
 wr(bind(boundary),'residual_size_boundary_cases.csv')
 stopifnot(all(vapply(c(0,1e-12,.1,1,1e6),function(a)a/max(a,sqrt(.Machine$double.eps))<=1,logical(1))))
 t<-paired(c(0,1,1,0),c(0,0,1,1),TRUE);stopifnot(t['mean']==0,t['n_up']==1,t['n_down']==1,t['lower']<0,t['upper']>0)
 t<-ci(c(NA,NA));stopifnot(t['n_finite']==0,is.na(t['mean']));t<-ci(c(0,NA,1));stopifnot(t['n_finite']==2,t['mean']==.5)
 wr(data.frame(check=c('235 immutable input SHA256 checks','98800 OK records;282 frozen cells','All data/algorithm/null seed formulas','94800 TSR records: realizedR0>=minPts; no infeasible reason','Exact TCLUST trimmed budget in every end-to-end TSR row','20 controlled contrast proofs including macro fit/count identity','Metric-specific finite masks and paired uncertainty fixtures','Analytic equality density<=1 arithmetic fixture'),status='PASS'),'validation_checks.csv','evidence')
 writeLines(capture.output(sessionInfo()),file.path(C,'evidence/R_session_info.txt'))
 # Plot only existing finite outcomes; no new data or model evaluation.
 sm<-bind(sums);pdf(file.path(C,'figures/minpts_existing_sensitivity.pdf'),width=9,height=6.2,useDingbats=FALSE)
 par(mfrow=c(2,2),mar=c(4,4,3,1));for(sc in c('positive','negative','weak','high_p25')){
  a<-sm[sm$scenario==sc&sm$density==2&sm$metric=='any_native_micro',];a<-a[order(a$minPts),]
  plot(a$minPts,a$mean,type='b',pch=19,col='#235985',ylim=c(0,1),xlab='minPts (max[5, multiplier x p])',ylab='Native Stage-II acceptance',main=paste(sc,'; p =',a$p[1]));arrows(a$minPts,a$lower,a$minPts,a$upper,angle=90,code=3,length=.04,col='#235985');abline(v=max(5,2*a$p[1]),lty=3,col='#666666');mtext('200 paired datasets per setting; exact pointwise 95% intervals',cex=.65)
 };dev.off()
 cat('COMMENT5_READONLY_AUDIT_OK;98800 records;20 contrasts; no model fits\n')
}
audit()
