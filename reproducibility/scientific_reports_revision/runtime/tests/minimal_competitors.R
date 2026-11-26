# A small deterministic API/fit smoke test, not a performance experiment.
source('R/load.R')
out<-Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','results/final_acceptance/local_20261119');dir.create(out,recursive=TRUE,showWarnings=FALSE)
set.seed(7015401)
x<-rbind(matrix(rnorm(80,sd=.18),40,2)+matrix(c(-3,0),40,2,byrow=TRUE),matrix(rnorm(80,sd=.18),40,2)+matrix(c(3,0),40,2,byrow=TRUE))
rows<-list()
for(method in c('OTRIMLE','tEIGEN')) {
 set.seed(7015402)
 fit<-if(method=='OTRIMLE')legacy$fit_otrimle_fixed(x,G=2) else legacy$fit_teigen_bic(x,Gs=2)
 stopifnot(length(fit$cluster)==nrow(x),all(is.finite(fit$cluster)),all(fit$cluster>=0))
 if(method=='OTRIMLE')stopifnot(is.null(fit$code)||fit$code!=0)
 rows[[method]]<-data.frame(method=method,n=nrow(x),groups=length(unique(fit$cluster[fit$cluster>0])),status='PASS')
}
write.csv(do.call(rbind,rows),file.path(out,'minimal_competitors.csv'),row.names=FALSE)
cat('MINIMAL_OTRIMLE_TEIGEN_OK\n')
