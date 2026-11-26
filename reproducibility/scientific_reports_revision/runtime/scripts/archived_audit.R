source('R/load.R')
files<-c('architecture','stress','negative_control');out<-list()
for(study in files) {
 z<-read.csv(gzfile(file.path(SOURCE_ROOT,'github_results/raw_results',paste0(study,'_replications.csv.gz'))),stringsAsFactors=FALSE)
 for(method in unique(z$method)) {
  w<-z[z$method==method,];if(study=='architecture')w<-w[w$scenario_code!='N0',]
  ok<-w$status=='OK'
  out[[length(out)+1]]<-data.frame(study=study,method=method,n_total=nrow(w),n_success=sum(ok),n_failed=sum(!ok),MajorARI=mean(w$MajorARI[ok],na.rm=TRUE),MinorityClusterF1=mean(w$MinorityClusterF1[ok],na.rm=TRUE),MinorityDetectionF1=mean(w$MinorityDetectionF1[ok],na.rm=TRUE),mean_QHat=mean(w$QHat[ok],na.rm=TRUE),QHat1=sum(w$QHat[ok]==1),MajorARI0=sum(w$MajorARI[ok]==0),requested_K=if(length(unique(w$requested_k))==1)w$requested_k[1] else NA)
 }
}
write.csv(bind_rows(out),'tables/archived_metrics_recomputed.csv',row.names=FALSE)
public<-c('helpers.R','tsrclust.R','DESCRIPTION')
cmp<-bind_rows(lapply(public,function(f) {
 local<-file.path(SOURCE_ROOT,'RPackageNew/TSRClust',if(f=='DESCRIPTION')f else file.path('R',f));pub<-file.path('config/public_snapshot',f)
 data.frame(file=f,local_md5=unname(tools::md5sum(local)),public_md5=unname(tools::md5sum(pub)),same_text=identical(readLines(local,warn=FALSE),readLines(pub,warn=FALSE)))
}))
write.csv(cmp,'reports/public_local_alignment.csv',row.names=FALSE);print(cmp)
