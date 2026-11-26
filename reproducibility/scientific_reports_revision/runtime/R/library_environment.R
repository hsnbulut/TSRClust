# Explicit per-process library initialization; no package is installed or modified.
initialize_revision_library <- function(path=Sys.getenv('TSR_R_LIB'),required=FALSE) {
 if(!nzchar(path)) {
  candidate<-Sys.getenv('R_LIBS_USER')
  if(nzchar(candidate)&&dir.exists(path.expand(candidate)))path<-candidate
 }
 if(!nzchar(path)) {
  if(required)stop('Set TSR_R_LIB to the existing isolated R 4.5.2 library')
  return(invisible(NULL))
 }
 path<-normalizePath(path.expand(path),mustWork=TRUE)
 if(!dir.exists(path))stop('Isolated library is not a directory: ',path)
 # R_LIBS_USER affects startup only; .libPaths must also be set in THIS process.
 Sys.setenv(TSR_R_LIB=path,R_LIBS_USER=path,R_LIBS=path,R_LIBS_SITE='NULL',
            R_ENVIRON='/dev/null',R_ENVIRON_USER='/dev/null',
            R_PROFILE='/dev/null',R_PROFILE_USER='/dev/null')
 .libPaths(c(path,.Library),include.site=FALSE)
 if(!identical(normalizePath(.libPaths()[1L]),path))stop('Could not initialize isolated library')
 invisible(path)
}

# Metadata preflight is cheap and read-only; namespace loading is verified separately.
assert_revision_environment <- function(lockfile='config/dependencies/source_lock.csv') {
 if(as.character(getRversion())!='4.5.2')stop('Actual R must be exactly 4.5.2')
 lib<-initialize_revision_library(required=TRUE)
 lock<-read.csv(lockfile,stringsAsFactors=FALSE)
 for(i in seq_len(nrow(lock))) {
  p<-lock$package[i];desc<-file.path(lib,p,'DESCRIPTION')
  if(!file.exists(desc))stop('Pinned package missing from isolated library: ',p)
  d<-read.dcf(desc)
  if(d[1,'Version']!=lock$version[i])stop('Pinned version mismatch: ',p)
  if('NeedsCompilation'%in%colnames(d)&&d[1,'NeedsCompilation']=='yes') {
   if(!'Built'%in%colnames(d)||!grepl('R 4[.]5[.]',d[1,'Built'])||!grepl(R.version$platform,d[1,'Built'],fixed=TRUE))stop('Native package R/architecture mismatch: ',p)
  }
 }
 invisible(TRUE)
}
