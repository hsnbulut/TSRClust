# Reviewer 1 Comment 2. Read immutable scalar results; no model fits/simulation.
stopifnot(as.character(getRversion())=='4.5.2')
f<-grep('^--file=',commandArgs(),value=TRUE)[1];code<-dirname(normalizePath(sub('^--file=','',f)))
c2<-dirname(code);rev<-dirname(dirname(c2));project<-dirname(rev)
write<-function(x,name,where='tables')write.csv(x,file.path(c2,where,name),row.names=FALSE,na='NA')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
inv<-read.csv(file.path(c2,'evidence/source_inventory.csv'))
stopifnot(all(vapply(seq_len(nrow(inv)),function(i)sha(inv$path[i])==inv$sha256[i],logical(1))))
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
stopifnot(sha(archive)=='c1852811d7dd94cf8f7b2dbae66d6c17f4284a6cdcf6abbfd8721d73076f94f6')
work<-tempfile('c2_archive_');dir.create(work)
run<-function(){
 on.exit(unlink(work,recursive=TRUE))
 untar(archive,files='replications.csv',exdir=work)
 x<-read.csv(file.path(work,'replications.csv'),check.names=FALSE)
 d<-read.csv(file.path(project,'revision_scientific_reports/config/design.csv'));d<-d[d$tier==1,]
 stopifnot(nrow(x)==98800,nrow(d)==282,all(x$status=='OK'),!anyDuplicated(x[c('cell_id','replication')]))
 for(n in names(d))stopifnot(identical(as.character(x[[n]]),as.character(d[[n]][match(x$cell_id,d$cell_id)])))
 x$active_novelty<-x$method=='TSRClust'&x$mode=='end_to_end'
 vals<-do.call(rbind,lapply(sort(unique(d$novelty)),function(v)data.frame(novelty=v,design_cells=sum(d$novelty==v),all_fit_records=sum(x$novelty==v),active_end_to_end_records=sum(x$novelty==v&x$active_novelty),residual_only_records=sum(x$novelty==v&x$mode=='residual_only'),other_method_records=sum(x$novelty==v&x$method!='TSRClust'))))
 write(vals,'evaluated_novelty_values.csv')
 quant<-expand.grid(p=sort(unique(d$p)),novelty=sort(unique(d$novelty)));quant$squared_cutoff<-qchisq(quant$novelty,quant$p);quant$radial_cutoff<-sqrt(quant$squared_cutoff)
 write(quant,'chi_square_cutoffs.csv')
 factors<-do.call(rbind,lapply(c('requested_alpha','structure_alpha','minority_proximity','minority_shape','separation','p','noise_type','reference'),function(n) {
  z<-as.data.frame(table(x[[n]]));names(z)<-c('value','fit_records');z$factor<-n;z[,c('factor','value','fit_records')]
 }));write(factors,'factor_coverage.csv')
 metrics<-c('AnyMicroDiscovery','MinorityClusterF1','MinorityDetectionF1_all_extra_v2','MajorARI','SpuriousExtraClusters')
 # A matching stratum includes EVERY frozen design field except threshold,
 # cell identifier and repetition count. No cross-reference/rule pooling.
 keys<-setdiff(names(d),c('cell_id','novelty','reps'))
 key<-apply(d[,keys],1,paste,collapse='|')
 groups<-split(seq_len(nrow(d)),key)
 groups<-Filter(function(ii)length(unique(d$novelty[ii]))>1&&any(d$novelty[ii]==.995),groups)
 stopifnot(length(groups)>0)
 scalar<-function(v,binary=FALSE) {
  v<-v[is.finite(v)];n<-length(v)
  if(!n)return(c(n_finite=0,mean=NA,mcse=NA,lo=NA,hi=NA))
  m<-mean(v);se<-if(n>1)sd(v)/sqrt(n) else NA_real_
  ci<-if(binary)binom.test(sum(v),n)$conf.int else if(n>1)m+c(-1,1)*qt(.975,n-1)*se else c(NA,NA)
  c(n_finite=n,mean=m,mcse=se,lo=ci[1],hi=ci[2])
 }
 paired<-list();groupmeta<-list();cellstats<-list();i<-0;j<-0;g<-0
 for(ii in groups) {
  g<-g+1;ds<-d[ii,];basecell<-ds$cell_id[ds$novelty==.995];stopifnot(length(basecell)==1)
  groupmeta[[g]]<-cbind(group_id=g,ds)
  b<-x[x$cell_id==basecell,]
  for(cid in ds$cell_id){
   a<-x[x$cell_id==cid,];dr<-d[d$cell_id==cid,]
   for(metric in metrics){
    j<-j+1;z<-scalar(a[[metric]],metric=='AnyMicroDiscovery')
    cellstats[[j]]<-cbind(dr,metric=metric,n_records=nrow(a),as.data.frame(as.list(z)),interval=if(metric=='AnyMicroDiscovery')'exact binomial 95%' else 'pointwise t Monte Carlo 95%')
   }
   if(cid==basecell)next
   m<-merge(a,b,by='replication',suffixes=c('_alt','_base'))
   stopifnot(nrow(m)==min(nrow(a),nrow(b)),!anyDuplicated(m$replication))
   for(seed in c('data_seed','algorithm_seed','null_seed'))stopifnot(all(m[[paste0(seed,'_alt')]]==m[[paste0(seed,'_base')]]))
   for(metric in metrics){
    va<-m[[paste0(metric,'_alt')]];vb<-m[[paste0(metric,'_base')]];ok<-is.finite(va)&is.finite(vb);delta<-va[ok]-vb[ok];z<-scalar(delta)
    plus<-minus<-NA_integer_;typ<-'pointwise paired t Monte Carlo 95%; zero variance gives degenerate empirical interval'
    if(metric=='AnyMicroDiscovery') {
     plus<-sum(delta==1);minus<-sum(delta== -1);n<-length(delta)
     cp<-binom.test(plus,n,conf.level=.975)$conf.int;cm<-binom.test(minus,n,conf.level=.975)$conf.int
     z[c('lo','hi')]<-c(cp[1]-cm[2],cp[2]-cm[1]);typ<-'conservative exact discordance bounds, Bonferroni 95%'
    }
    i<-i+1;paired[[i]]<-cbind(group_id=g,alternative_cell=cid,baseline_cell=basecell,experiment=dr$experiment,scenario=dr$scenario,reference=dr$reference,rule=dr$rule,novelty=dr$novelty,baseline_novelty=.995,p=dr$p,requested_alpha=dr$requested_alpha,minority_proximity=dr$minority_proximity,minority_ratio=dr$minority_ratio,noise_type=dr$noise_type,metric=metric,n_joined=nrow(m),as.data.frame(as.list(z)),positive_discordant=plus,negative_discordant=minus,interval=typ)
   }
  }
 }
 pair<-do.call(rbind,paired);cells<-do.call(rbind,cellstats);gm<-do.call(rbind,groupmeta)
 write(pair,'paired_novelty_differences.csv');write(cells,'sensitivity_cell_performance.csv');write(gm,'paired_design_strata.csv')
 scientific_keys<-setdiff(names(d),c('cell_id','experiment','scenario','reps'))
 equivalence<-apply(gm[,scientific_keys],1,paste,collapse='|')
 aliases<-gm[duplicated(equivalence)|duplicated(equivalence,fromLast=TRUE),c('group_id','cell_id','experiment','scenario','novelty','dataset_id')]
 write(aliases,'duplicate_scientific_configurations.csv','evidence')
 # Context is cell-specific; differences in other factors are not threshold effects.
 context<-x[x$active_novelty&x$novelty==.995,]
 cr<-list();j<-0
 for(cid in unique(context$cell_id))for(metric in metrics){
  a<-context[context$cell_id==cid,];dr<-d[d$cell_id==cid,];j<-j+1
  cr[[j]]<-cbind(dr,metric=metric,n_records=nrow(a),as.data.frame(as.list(scalar(a[[metric]],metric=='AnyMicroDiscovery'))))
 }
 write(do.call(rbind,cr),'default_context_performance.csv')
 # Compact manuscript subset: actual baseline-reference B proximity cells, separate sizes.
 manuscript<-cells[cells$experiment=='B'&cells$reference=='legacy_empirical'&cells$novelty%in%c(.95,.995,.999)&cells$metric%in%c('AnyMicroDiscovery','MinorityClusterF1','MinorityDetectionF1_all_extra_v2'),]
 write(manuscript,'proximity_manuscript_evidence.csv')
 # Proposed gap-filling only. Original Tier1 manifest is never edited.
 gap<-d[d$experiment=='A-END'&d$rule=='legacy_rule'&d$reference=='legacy_empirical'&d$p==5&d$requested_alpha==.1&d$noise_type%in%c('correlated_gaussian','heterogeneous'),]
 stopifnot(nrow(gap)==2)
 proposal<-do.call(rbind,lapply(seq_len(nrow(gap)),function(h)do.call(rbind,lapply(c(.95,.975,.99,.999),function(v)cbind(proposal_id=paste0('C2-gap-',gap$cell_id[h],'-',v),baseline_cell=gap$cell_id[h],proposed_novelty=v,rep_start=1,rep_end=100,additional_fits=100,gap[h,setdiff(names(gap),c('cell_id','novelty','reps'))])))))
 write(proposal,'optional_gap_proposal_NOT_EXECUTED.csv')
 overview<-data.frame(item=c('successful_records','design_cells','evaluated_thresholds','paired_recorded_strata','paired_unique_scientific_strata','paired_alternative_cells','paired_metric_rows','paired_joined_fits_not_independent','proposed_additional_fits'),value=c(nrow(x),nrow(d),nrow(vals),length(groups),length(unique(equivalence))/length(unique(gm$novelty)),length(unique(pair$alternative_cell)),nrow(pair),sum(pair$n_joined[pair$metric=='AnyMicroDiscovery']),sum(proposal$additional_fits)))
 write(overview,'audit_counts.csv')
 saveRDS(list(counts=overview,values=vals,paired=pair,cells=cells,context=do.call(rbind,cr)),file.path(c2,'evidence/compact_novelty_audit.rds'))
 # Exact implementation checks on fixed matrices, not empirical simulations.
 e<-new.env();sys.source(file.path(project,'Scientific Reports Manuscript/reproducibility/core/core_tsr.R'),e)
 a<-cbind(c(-2,-1,1,2,5,6,8,9),c(-4,-2,2,4,1,2,4,5));lab<-c(rep(1,4),rep(2,4));ref<-e$make_macro_reference(a,lab,2)
 S<-cov(a[lab==1,,drop=FALSE]);ridge<-max(mean(diag(S)),1)*1e-8
 stopifnot(max(abs(ref[[1]]$covariance-S-diag(ridge,2)))<1e-12,all(eigen(ref[[1]]$covariance,symmetric=TRUE)$values>0))
 candidate<-rbind(c(20,30),c(22,32));manual<-min(vapply(ref,function(r)mahalanobis(colMeans(candidate),r$center,r$covariance),numeric(1)))
 stopifnot(abs(e$candidate_novelty_distance(candidate,1:2,ref)-manual)<1e-8,is.infinite(e$candidate_novelty_distance(candidate,1:2,NULL)))
 fallback<-e$make_macro_reference(a,rep(1,nrow(a)),2)[[2]];stopifnot(identical(fallback$center,colMeans(a)))
 checks<-data.frame(check=c('immutable sources and archive','complete 98800 frozen records','all design columns agree','all threshold pairs match data algorithm and null seeds','legacy empirical reference ridge and centroid minimum','empty macro reference returns Inf and empty group fallback'),status='PASS',detail=c('75 source hashes and archive SHA256 checked','282 cells; unique cell/replication; status OK','No frozen parameter changed',paste(length(unique(pair$alternative_cell)),'alternative cells;',sum(pair$n_joined[pair$metric=='AnyMicroDiscovery']),'joined pairs'),'Fixed deterministic matrices; no model fitting','Residual-only novelty is inactive; fallback uses whole analysis mean/identity'))
 write(checks,'validation_checks.csv','evidence');writeLines(capture.output(sessionInfo()),file.path(c2,'evidence/R_session_info.txt'))
 print(overview);cat('COMMENT2_READONLY_AUDIT_OK\n')
}
run()
