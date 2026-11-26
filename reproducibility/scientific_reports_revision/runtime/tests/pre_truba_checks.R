source('R/load.R')
checks<-list();test<-function(name,expr) {force(expr);checks[[length(checks)+1L]]<<-data.frame(check=name,status='PASS');cat('PASS',name,'\n')}
test('Actual R 4.5.2 and all required APIs load',{
 stopifnot(as.character(getRversion())=='4.5.2',requireNamespace('otrimle'),requireNamespace('teigen'),requireNamespace('robustbase'))
})
test('same observed/null grids for legacy and p-only; recalibration regenerates',{
 original<-legacy$make_eps_grid;calls<-list()
 legacy$make_eps_grid<-function(x,min_pts,grid_size=31L) {
  ans<-original(x,min_pts,grid_size);calls[[length(calls)+1L]]<<-ans;ans
 }
 tryCatch({
  set.seed(123);xx<-rbind(matrix(rnorm(40,sd=.03),20,2),matrix(runif(160,-4,4),80,2))
  calls<-list();set.seed(42);a<-revision_micro(xx,5,B=3L,alpha=.5,density=1.01);ca<-calls
  calls<-list();set.seed(42);b<-revision_micro(xx,5,B=3L,alpha=.5,density=1.01,rule='pvalue_only');cb<-calls
  calls<-list();set.seed(42);c<-revision_micro(xx,5,B=3L,alpha=.5,density=1.01,rule='recalibrated_rule',conditional_allowed=TRUE);cc<-calls
  stopifnot(length(ca)==1,length(cb)==1,length(cc)==4,identical(a$null_scores,b$null_scores),identical(ca,cb),any(vapply(cc[-1],function(z)!identical(z,cc[[1]]),logical(1))))
 },finally=legacy$make_eps_grid<-original)
})
test('legacy acceptance remains submitted residual rule',{
 set.seed(77);xx<-rbind(matrix(rnorm(40,sd=.03),20,2),matrix(runif(160,-4,4),80,2));rng<-.Random.seed
 a<-legacy$discover_validated_microclusters(xx,5,49,.05,2,NULL,qchisq(.995,2));.Random.seed<-rng
 b<-revision_micro(xx,5,49,ref=NULL)
 stopifnot(identical(a$cluster,b$cluster),isTRUE(all.equal(a$p_value,b$p_value)),isTRUE(all.equal(a$threshold,b$threshold)))
})
test('oracle trimming uses pure-noise proportion once',{
 design<-read.csv('config/design.csv');row<-design[design$method=='Oracle_TCLUST'&design$scenario=='positive',][1,]
 z<-fit_revision(row,1);stopifnot(z$row$status=='OK')
 stopifnot(row$requested_alpha==row$residual_total-row$minority_ratio,
  sum(z$diagnostics$raw_cluster==0)==floor(row$n*row$requested_alpha))
})
test('explicit metric aliases and macro contaminants',{
 y<-c(rep(1,100),rep(2,100),rep(3,20));pred<-c(rep(1,90),rep(3,10),rep(2,100),rep(3,20))
 z<-revision_metrics(y,pred,2)$metrics
 stopifnot(z$MinorityClusterF1_legacy_v1==.8,abs(z$MinorityPrecision_legacy_v1-2/3)<1e-12,
  abs(z$MinorityPrecision_all_extra_v2-2/3)<1e-12,z$MinorityDetectionF1_all_extra_v2==.8)
})
test('DP agrees with exhaustive optional assignment',{
 brute<-function(M,i=1,used=integer()) {
  if(i>nrow(M))return(0)
  max(c(brute(M,i+1,used),vapply(setdiff(seq_len(ncol(M)),used),function(j)M[i,j]+brute(M,i+1,c(used,j)),numeric(1))))
 }
 for(seed in 1:10){set.seed(seed);M<-matrix(sample(c(0,.2,.5,1),12,TRUE),3,4);path<-legacy$max_assignment_dp(M);score<-sum(vapply(seq_len(3),function(i)if(path[i])M[i,path[i]] else 0,numeric(1)));stopifnot(abs(score-brute(M))<1e-12)}
})
write.csv(do.call(rbind,checks),file.path(Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','results/pre_truba'),'targeted_checks.csv'),row.names=FALSE)
cat('PRE_TRUBA_TARGETED_CHECKS_OK',length(checks),'\n')
