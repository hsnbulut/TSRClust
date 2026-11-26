# Record every exported field and its finite denominator for OTRIMLE.
f<-grep('^--file=',commandArgs(),value=TRUE)[1];c3<-dirname(dirname(normalizePath(sub('^--file=','',f))));project<-dirname(dirname(dirname(c3)))
w<-tempfile('c3_columns_');dir.create(w)
column_audit<-function(){
 on.exit(unlink(w,recursive=TRUE));out<-list()
 for(st in c('architecture','stress','negative_control','tier1')){
  if(st=='tier1'){
   untar(Sys.getenv('TSR_TIER1_ARCHIVE'),files='replications.csv',exdir=w)
   x<-read.csv(file.path(w,'replications.csv'))
  }else x<-readRDS(file.path(project,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined',paste0(st,'_replications.rds')))
  x<-x[x$method=='OTRIMLE_fixedK',]
  for(nm in names(x))out[[length(out)+1]]<-data.frame(source=st,field=nm,type=class(x[[nm]])[1],n_records=nrow(x),n_missing=sum(is.na(x[[nm]])),n_finite=if(is.numeric(x[[nm]]))sum(is.finite(x[[nm]])) else NA_integer_)
 }
 write.csv(do.call(rbind,out),file.path(c3,'evidence/export_field_availability.csv'),row.names=FALSE)
 cat('ALL_EXPORT_FIELD_FINITE_DENOMINATORS_RECORDED\n')
}
column_audit()
