# Structural deduction from stored component budget and matching counts; no fits.
f<-grep('^--file=',commandArgs(),value=TRUE)[1];c3<-dirname(dirname(normalizePath(sub('^--file=','',f))));project<-dirname(dirname(dirname(c3)))
w<-tempfile('c3_capacity_');dir.create(w)
main<-function(){
 on.exit(unlink(w,recursive=TRUE));out<-list();j<-0
 for(st in c('architecture','stress')){
  x<-readRDS(file.path(project,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined',paste0(st,'_replications.rds')))
  x<-x[x$method=='OTRIMLE_fixedK'&x$q>0,];x$G<-x$selected_g
  case<-is.finite(x$G)&x$G==2&x$QHat==1&x$MajorRetention>0
  z<-x[case,];stopifnot(all(z$main_k>=2));j<-j+1
  out[[j]]<-data.frame(source='Original submitted',cohort=st,n_positive_fits=nrow(x),n_single_macro_component_identified=nrow(z),fraction_identified=nrow(z)/nrow(x),mean_macro_fraction_in_that_component=mean(z$MajorRetention),mean_MajorARI_identified=mean(z$MajorARI),mean_MinorityClusterF1_identified=mean(z$MinorityClusterF1),basis='Positive labels <= G=2; one eligible extra; positive macro retention implies exactly one macro-matched positive label. No membership or confusion count reconstruction.')
 }
 untar(Sys.getenv('TSR_TIER1_ARCHIVE'),files='replications.csv',exdir=w);x<-read.csv(file.path(w,'replications.csv'));x<-x[x$method=='OTRIMLE_fixedK'&x$q>0,]
 for(s in unique(x$scenario)){
  a<-x[x$scenario==s,];case<-a$estimated_K==2&a$QHat==1&a$MajorRetention>0;z<-a[case,];stopifnot(all(z$main_k>=2));j<-j+1
  out[[j]]<-data.frame(source='Tier1 revision',cohort=s,n_positive_fits=nrow(a),n_single_macro_component_identified=nrow(z),fraction_identified=nrow(z)/nrow(a),mean_macro_fraction_in_that_component=if(nrow(z))mean(z$MajorRetention) else NA_real_,mean_MajorARI_identified=if(nrow(z))mean(z$MajorARI) else NA_real_,mean_MinorityClusterF1_identified=if(nrow(z))mean(z$MinorityClusterF1) else NA_real_,basis='Positive labels <= G=2; one eligible extra; positive macro retention implies exactly one macro-matched positive label. No membership or confusion count reconstruction.')
 }
 write.csv(do.call(rbind,out),file.path(c3,'tables/component_capacity_diagnostics.csv'),row.names=FALSE);print(do.call(rbind,out))
}
main()
