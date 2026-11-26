# Deterministic STATISTIC fixtures, not generated datasets or new scientific fits.
a<-grep('^--file=',commandArgs(),value=TRUE)[1];C<-dirname(dirname(normalizePath(sub('^--file=','',a))));R<-dirname(dirname(C));P<-dirname(R)
stopifnot(as.character(getRversion())=='4.5.2')
module<-new.env();sys.source(file.path(P,'revision_scientific_reports/R/validation.R'),module)
orig<-new.env();sys.source(file.path(P,'Scientific Reports Manuscript/reproducibility/core/core_tsr.R'),orig)
pkg<-new.env();sys.source(file.path(P,'Scientific Reports Manuscript/RPackageNew/TSRClust/R/helpers.R'),pkg)
# Independent reference compares integers, without calling production mc_p/mc_accept.
reference<-function(s,n,a,strict){r<-sum(n>=s);isTRUE(s>0&&(r+1)<=a*(length(n)+1)&&(!strict||all(s>n)))}
truth<-list();exhaustive<-list();k<-0L
for(B in c(1L,9L,18L,19L,49L,99L,199L))for(alpha in c(0,.005,.01,.019,.02,.04,.05,.06,.1,1))for(e in 0:B){
 n<-c(rep(11,e),rep(1,B-e));s<-10
 l<-reference(s,n,alpha,TRUE);p<-reference(s,n,alpha,FALSE)
 stopifnot(identical(l,module$mc_accept(s,n,alpha,'legacy_rule')),identical(p,module$mc_accept(s,n,alpha,'pvalue_only')))
 k<-k+1L;exhaustive[[k]]<-data.frame(B=B,structure_alpha=alpha,exceedances_ge=e,p_MC=(1+e)/(B+1),pvalue_only=p,legacy=l)
}
# Execute the unmodified full ORIGINAL and PACKAGE decision bodies. Collaborators
# returning grids/candidates/null matrices are replaced only in disposable environments.
# This isolates decision logic without dbscan/TCLUST fits or random observations.
actual<-function(e,prefix,obs,null,alpha=.05,n=6L,empty=FALSE,short_grid=FALSE){
 f<-get(paste0(prefix,'discover_validated_microclusters'),e)
 e2<-new.env(parent=e);environment(f)<-e2;calls<-0L
 cand<-function(scores){if(!length(scores))return(list());list(list(eps=1,labels=as.integer(c(rep(seq_along(scores),each=3),rep(0L,max(0L,n-3L*length(scores))))),scores=data.frame(cluster=seq_along(scores),size=rep(3L,length(scores)),density_ratio=rep(2,length(scores)),score=scores)))}
 assign(paste0(prefix,'make_eps_grid'),function(...)if(short_grid)numeric() else 1,e2)
 assign(paste0(prefix,'simulate_diffuse_noise'),function(x)x,e2)
 assign(paste0(prefix,'scan_dbscan_candidates'),function(...){calls<<-calls+1L;if(short_grid)return(list());if(calls==1L)return(if(empty)list() else cand(obs));s<-null[calls-1L];if(is.na(s))return(cand(s));if(s<=0)return(list());cand(s)},e2)
 tryCatch({z<-f(matrix(0,n,2),3L,length(null),alpha,2,NULL,0);stopifnot(length(z$cluster)==n);list(status='OK',accepted=any(z$cluster>0),p=z$p_value,nlabels=sum(z$cluster>0),calls=calls,message='')},error=function(err)list(status='ERROR',accepted=NA,p=NA,nlabels=NA,calls=calls,message=conditionMessage(err)))
}
fixtures<-list(
 list('B49_zero_exceedances',49,.05,10,rep(1,49)),
 list('B49_one_greater',49,.05,10,c(11,rep(1,48))),
 list('B49_two_greater',49,.05,10,c(11,11,rep(1,47))),
 list('B49_one_exact_tie',49,.05,10,c(10,rep(1,48))),
 list('B49_two_exact_ties',49,.05,10,c(10,10,rep(1,47))),
 list('B49_all_ties',49,.05,10,rep(10,49)),
 list('B49_zero_observed_and_null',49,.05,0,rep(0,49)),
 list('B49_zero_observed_positive_null',49,.05,0,rep(1,49)),
 list('B49_positive_observed_zero_null',49,.05,10,rep(0,49)),
 list('B49_alpha_below_minimum',49,.019,10,rep(1,49)),
 list('B49_alpha_at_minimum',49,.02,10,rep(1,49)),
 list('B49_alpha04_one_exceedance',49,.04,10,c(11,rep(1,48))),
 list('B19_default_at_minimum',19,.05,10,rep(1,19)),
 list('B18_minimum_too_large',18,.05,10,rep(1,18)),
 list('B99_four_exceedances',99,.05,10,c(rep(11,4),rep(1,95))),
 list('B99_zero_exceedances',99,.05,10,rep(1,99)),
 list('B199_nine_exceedances',199,.05,10,c(rep(11,9),rep(1,190))),
 list('B199_ten_exceedances',199,.05,10,c(rep(11,10),rep(1,189))),
 list('B199_zero_exceedances',199,.05,10,rep(1,199)))
for(f in fixtures){
 id<-f[[1]];B<-f[[2]];alpha<-f[[3]];s<-f[[4]];null<-f[[5]];r<-sum(null>=s);l<-reference(s,null,alpha,TRUE);p<-reference(s,null,alpha,FALSE)
 oa<-actual(orig,'',s,null,alpha);pa<-actual(pkg,'.',s,null,alpha)
 stopifnot(oa$status=='OK',pa$status=='OK',identical(oa$accepted,l),identical(pa$accepted,l),module$mc_accept(s,null,alpha,'pvalue_only')==p)
 if(s>0)stopifnot(oa$p==(1+r)/(B+1),pa$p==(1+r)/(B+1)) else stopifnot(oa$p==1,pa$p==1,oa$calls==1,pa$calls==1)
 truth[[length(truth)+1L]]<-data.frame(case=id,B=B,structure_alpha=alpha,observed_score=s,exceedances_ge=r,ties=sum(null==s),p_MC_formula=(1+r)/(B+1),pvalue_only=p,legacy=l,original_body=oa$accepted,package_body=pa$accepted,stored_global_p=oa$p,null_scans_executed=max(0,oa$calls-1L))
}
write.csv(do.call(rbind,truth),file.path(C,'tables/decision_truth_table.csv'),row.names=FALSE)
write.csv(do.call(rbind,exhaustive),file.path(C,'tables/exhaustive_boundary_checks.csv'),row.names=FALSE)
# A global pass does not retain every candidate: score10 passes; score2 tied with
# max(null)=2 fails strict candidate filtering. The unmodified body must keep3/6.
edge<-list()
for(kind in c('original','package')){
 en<-if(kind=='original')orig else pkg;pr<-if(kind=='original')'' else '.'
 z<-actual(en,pr,c(10,2),rep(2,49));stopifnot(z$status=='OK',z$p==.02,z$nlabels==3)
 edge[[length(edge)+1]]<-data.frame(implementation=kind,case='global_pass_lower_candidate_tied',status=z$status,accepted=z$accepted,stored_p=z$p,null_scans=max(0,z$calls-1),retained_observations=z$nlabels,message=z$message)
 for(case in c('empty_scan','short_residual','empty_grid','missing_observed','missing_null','infinite_observed')){
  ss<-if(case=='missing_observed')NA_real_ else if(case=='infinite_observed')Inf else 10
  ns<-if(case=='missing_null')c(NA_real_,rep(1,48)) else rep(1,49)
  z<-actual(en,pr,ss,ns,n=if(case=='short_residual')2L else 6L,empty=case=='empty_scan',short_grid=case=='empty_grid')
  if(case%in%c('empty_scan','short_residual','empty_grid'))stopifnot(z$status=='OK',!z$accepted,z$p==1)
  if(grepl('missing',case))stopifnot(z$status=='ERROR')
  if(case=='infinite_observed')stopifnot(z$status=='OK',z$accepted) # no is.finite guard in this decision body
  edge[[length(edge)+1]]<-data.frame(implementation=kind,case=case,status=z$status,accepted=z$accepted,stored_p=z$p,null_scans=max(0,z$calls-1),retained_observations=z$nlabels,message=z$message)
 }
}
write.csv(do.call(rbind,edge),file.path(C,'tables/edge_and_invalid_input_checks.csv'),row.names=FALSE,na='NA')
# General attainable boundary for nonnegative finite scores, B>=1, 0<=a<=1.
general<-expand.grid(B=c(9,18,19,49,99,199),structure_alpha=c(.005,.01,.019,.02,.04,.05,.1),KEEP.OUT.ATTRS=FALSE)
general$minimum_p<-1/(general$B+1);general$pvalue_max_exceedances<-mapply(function(B,a){ok<-(1+0:B)/(B+1)<=a;if(!any(ok))NA_integer_ else max((0:B)[ok])},general$B,general$structure_alpha)
general$strict_gate_can_accept<-general$minimum_p<=general$structure_alpha;general$legacy_accepted_p<-ifelse(general$strict_gate_can_accept,general$minimum_p,NA_real_)
general$ideal_continuous_rank_probability_strict<-ifelse(general$strict_gate_can_accept,general$minimum_p,0)
general$ideal_rank_bound_pvalue<-ifelse(is.na(general$pvalue_max_exceedances),0,(general$pvalue_max_exceedances+1)/(general$B+1))
write.csv(general,file.path(C,'tables/general_B_alpha_boundaries.csv'),row.names=FALSE,na='NA')
writeLines(c('PASS',paste('Exhaustive scalar cases:',k),paste('Original/package full-body fixtures:',length(fixtures)*2),paste('Edge/invalid-input full-body checks:',nrow(do.call(rbind,edge)))) ,file.path(C,'evidence/boundary_test_receipt.txt'))
cat('MC_BOUNDARY_TEST_PASS;',k,'exhaustive cases;',length(fixtures)*2,'full-body valid fixtures;',length(edge),'edge fixtures\n')
