source('R/load.R')
checks<-list();test<-function(name,expr) {force(expr);checks[[length(checks)+1L]]<<-data.frame(test=name,status='PASS');cat('PASS',name,'\n')}
error<-function(expr)inherits(tryCatch({force(expr);NULL},error=function(e)e),'error')
test('discrete p grid and maximum gate',{
 null<-1:49;stopifnot(mc_p(48.5,null)==.04,!mc_accept(48.5,null),mc_accept(48.5,null,rule='pvalue_only'),mc_p(50,null)==.02,mc_accept(50,null),!mc_accept(0,rep(0,49)))
})
test('rank rejection exact enumeration',{
 for(B in c(49,99,199)) {
  scores<-seq_len(B+1);legacy_rej<-p_rej<-0
  for(i in seq_along(scores)){ns<-scores[-i];legacy_rej<-legacy_rej+mc_accept(scores[i],ns);p_rej<-p_rej+mc_accept(scores[i],ns,rule='pvalue_only')}
  stopifnot(legacy_rej==1,p_rej==floor(.05*(B+1)))
 }
})
test('conditional samples preserve marginal extrema',{
 set.seed(1);x<-matrix(runif(100),50,2);y<-conditional_box(x)
 stopifnot(identical(apply(x,2,range),apply(y,2,range)),!identical(x,y),error(conditional_box(matrix(1,10,2))),error(revision_micro(x,5,rule='recalibrated_rule')))
})
test('nonoverlapping seed namespaces with paired configurations',{
 seeds<-unlist(lapply(0:2,function(stream)unlist(lapply(1:2,function(rr) vapply(1:2000,function(rep)seed_id(7,rep,stream,rr),integer(1))))))
 stopifnot(!anyDuplicated(seeds),seed_id(7,1,0)==seed_id(7,1,0),error(seed_id(1000,1,0)))
})
test('merged macros and genuine minority are compatible',{
 truth<-c(rep(1,100),rep(2,100),rep(3,20));pred<-c(rep(1,200),rep(2,20));z<-revision_metrics(truth,pred,2)$metrics
 stopifnot(z$MajorARI==0,z$MinorityClusterF1==1,z$MinorityPrecision==1,z$QHat==1)
})
test('macro fragments omitted by legacy detection denominator',{
 truth<-c(rep(1,100),rep(2,100),rep(3,20));pred<-c(rep(1,60),rep(4,40),rep(2,100),rep(3,20));z<-revision_metrics(truth,pred,2)$metrics
 stopifnot(z$MinorityDetectionF1==1,z$MinorityDetectionF1_all_extra_v2==.5,z$MinorityPrecision_all_extra_v2==1/3,z$MacroFalsePositiveReassignment==.2)
 write.csv(z,file.path(Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','tables'),'deterministic_fragmentation_example.csv'),row.names=FALSE)
})
test('macro contaminants count against minority precision',{
 truth<-c(rep(1,100),rep(2,100),rep(3,20));pred<-c(rep(1,90),rep(3,10),rep(2,100),rep(3,20));z<-revision_metrics(truth,pred,2)$metrics
 stopifnot(abs(z$MinorityPrecision-2/3)<1e-12,abs(z$MinorityClusterF1-.8)<1e-12)
})
test('metric invariance and absent minority',{
 y<-c(rep(1,20),rep(2,20),rep(3,10));p<-c(rep(8,20),rep(9,20),rep(7,10));a<-revision_metrics(y,p,2)$metrics;b<-revision_metrics(y,p*10,2)$metrics
 stopifnot(isTRUE(all.equal(a,b)),is.na(revision_metrics(rep(1:2,each=20),rep(1:2,each=20),2)$metrics$MinorityClusterF1),error(revision_metrics(1:2,1,2)))
})
test('DBSCAN feasibility follows chosen minPts',{
 x<-matrix(seq_len(35)*1e-4,7,5);stopifnot(all(dbscan::dbscan(x,eps=1,minPts=5)$cluster==1),revision_micro(x,10)$reason=='infeasible_minPts',length(legacy$make_eps_grid(x,5))>0)
})
# Match aligned package, submitted analysis, and revision on fixed (not searched) seeds.
pkg<-new.env(parent=globalenv());for(f in c('helpers.R','simulate.R','tsrclust.R'))source(file.path(SOURCE_ROOT,'RPackageNew/TSRClust/R',f),local=pkg)
aligned<-list()
for(sc in c('N0','M1','M3','M4')) {
 d<-pkg$simulate_tsr_data(n=500,p=2,scenario=sc,structure='elliptical',seed=1)
 x<-legacy$robust_standardize(d$x);set.seed(1001)
 macro<-legacy$fit_tclust_fixed(x,2,.1,20,50,legacy$TSR_CONFIG)
 # Keep the pre-null RNG state exactly as in submitted package/reference alignment.
 rng<-.Random.seed
 old<-legacy$fit_tsr_from_macro(x,macro$cluster,2,TRUE,49)
 .Random.seed<-rng
 new<-revision_micro(x[old$residual_index,,drop=FALSE],5,49,ref=legacy$make_macro_reference(x,macro$cluster,2))
 set.seed(1001)
 pp<-pkg$tsrclust(d$x,k=2,alpha=.1,restr_fact=20,null_reps=49,nstart=50,niter1=3,niter2=20,nkeep=5,opt='HARD',equal.weights=FALSE,parallel=FALSE,store_x=FALSE,trace=0)
 test(paste('baseline package/analysis/revision',sc),{
  stopifnot(identical(old$micro$cluster,new$cluster),identical(as.integer(old$cluster),pp$cluster))
  for(nm in c('eps','p_value','score','threshold'))stopifnot(isTRUE(all.equal(old$micro[[nm]],new[[nm]])))
  if(sc!='N0')stopifnot(any(new$cluster>0))
 })
 aligned[[sc]]<-data.frame(scenario=sc,n_micro=length(unique(new$cluster[new$cluster>0])),p=new$p_value,score=new$score)
}
write.csv(bind_rows(aligned),file.path(Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','tables'),'local_baseline_alignment.csv'),row.names=FALSE)
test('direct TCLUST reference dimensions',{
 r<-revision_reference(x,macro,'tclust_direct');stopifnot(identical(r[[1]]$center,macro$centers[,1]),identical(r[[1]]$covariance,macro$cov[,,1]))
})
test('all array IDs map once to existing cells and rep blocks',{
 d<-read.csv('config/design.csv')
 for(f in c('pilot','tier1','tier2')) {
  t<-read.csv(paste0('config/tasks_',f,'.csv'));stopifnot(identical(t$task_id,seq_len(nrow(t))),all(t$cell_id%in%d$cell_id),all(t$rep_end>=t$rep_start))
  keys<-unlist(lapply(seq_len(nrow(t)),function(i)paste(t$cell_id[i],seq(t$rep_start[i],t$rep_end[i]),sep='_')));stopifnot(!anyDuplicated(keys))
 }
})
test('data and fit deterministic across repeated calls',{
 d<-read.csv('config/design.csv');r<-d[d$experiment=='C'&d$scenario=='positive',][1,]
 a<-fit_revision(r,1);b<-fit_revision(r,1);stopifnot(a$row$status=='OK',identical(a$diagnostics$predicted,b$diagnostics$predicted),identical(a$diagnostics$micro$null_scores,b$diagnostics$micro$null_scores))
 r$method<-'NOT_A_METHOD';z<-fit_revision(r,1);stopifnot(z$row$status=='ERROR',nzchar(z$row$error))
})
test('atomic result round trip',{
 p<-file.path(tempdir(),'tsr_roundtrip.rds');atomic_rds(list(x=1),p);stopifnot(readRDS(p)$x==1);unlink(p)
})
write.csv(bind_rows(checks),file.path(Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','reports'),'local_tests.csv'),row.names=FALSE)
cat('ALL_TESTS_OK',length(checks),'\n')
