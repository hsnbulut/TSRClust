# Fresh-process probe. Called by the regression test, never installs packages.
a<-commandArgs(TRUE);stopifnot(length(a)==2L)
mode<-a[1];outfile<-a[2];expected<-normalizePath(Sys.getenv('TSR_R_LIB'),mustWork=TRUE)
if(mode=='redirect') {
 empty<-tempfile('empty-library-');dir.create(empty)
 Sys.setenv(R_LIBS_USER=empty)
 .libPaths(c(empty,.Library),include.site=FALSE)
}
before<-list(R_LIBS_USER=Sys.getenv('R_LIBS_USER'),libPaths=.libPaths())
if(mode=='raw') {
 saveRDS(list(before=before,expected_visible=expected%in%normalizePath(.libPaths())),outfile)
 quit(save='no')
}
source('R/library_environment.R');lib<-initialize_revision_library(required=TRUE)
stopifnot(as.character(getRversion())=='4.5.2')
lock<-read.csv('config/dependencies/source_lock.csv',stringsAsFactors=FALSE)
checks<-lapply(seq_len(nrow(lock)),function(i) {
 p<-lock$package[i];ok<-requireNamespace(p,quietly=TRUE)
 path<-if(ok)getNamespaceInfo(asNamespace(p),'path') else ''
 desc<-if(ok)utils::packageDescription(p) else NULL
 data.frame(package=p,loaded=ok,version_matches=ok&&desc$Version==lock$version[i],
  isolated=ok&&normalizePath(dirname(path))==lib,library=path)
})
z<-do.call(rbind,checks)
saveRDS(list(before=before,after=list(R_LIBS_USER=Sys.getenv('R_LIBS_USER'),libPaths=.libPaths()),checks=z),outfile)
if(!all(z$loaded&z$version_matches&z$isolated))stop('Fresh child cannot see all pinned packages in isolated library')
cat('FRESH_CHILD_ALL_17_OK mode=',mode,'\n',sep='')
