a<-commandArgs(TRUE);if(!length(a))stop('Usage: verify_dependencies.R isolated_library')
if(as.character(getRversion())!='4.5.2')stop('Actual R must be 4.5.2')
source('R/library_environment.R')
lib<-initialize_revision_library(a[1],required=TRUE)
assert_revision_environment()
cat('R_LIBS_USER=',Sys.getenv('R_LIBS_USER'),'\n.libPaths=',paste(.libPaths(),collapse=' | '),'\n',sep='')
lock<-read.csv('config/dependencies/source_lock.csv',stringsAsFactors=FALSE,check.names=FALSE)
rows<-lapply(seq_len(nrow(lock)),function(i) {
 p<-lock$package[i];ok<-requireNamespace(p,quietly=TRUE);desc<-if(ok)packageDescription(p) else NULL
 path<-if(ok)getNamespaceInfo(asNamespace(p),'path') else ''
 data.frame(package=p,expected=lock$version[i],actual=if(ok)desc$Version else NA,loaded=ok,from_isolated_library=ok&&normalizePath(dirname(path))==lib,version_matches=ok&&desc$Version==lock$version[i],Built=if(ok)desc$Built else NA)
})
z<-do.call(rbind,rows);outdir<-Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','results/pre_truba');dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
write.csv(z,file.path(outdir,'dependency_verification.csv'),row.names=FALSE)
writeLines(capture.output(sessionInfo()),file.path(outdir,'session_R452.txt'))
print(z[c('package','actual','loaded','from_isolated_library','version_matches')])
stopifnot(all(z$loaded),all(z$from_isolated_library),all(z$version_matches))
# Explicit API arguments required by submitted wrappers.
stopifnot(all(c('data','G','erc','npr.max','ncores','monitor')%in%names(formals(otrimle::otrimle))),
 all(c('x','Gs','models','init','scale','verbose','parallel.cores')%in%names(formals(teigen::teigen))),
 all(c('x')%in%names(formals(robustbase::covMcd))))
cat('PINNED_DEPENDENCY_API_OK\n')
