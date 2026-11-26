# Linux installation must run on a node belonging to an active SLURM allocation.
if(identical(unname(Sys.info()[['sysname']]),'Linux')) {
 job<-Sys.getenv('SLURM_JOB_ID');nodes<-Sys.getenv('SLURM_JOB_NODELIST',Sys.getenv('SLURM_NODELIST'))
 if(!nzchar(job)||!nzchar(nodes)||!nzchar(Sys.which('scontrol')))stop('Linux package installation requires a SLURM compute-node allocation; do not install on the login node')
 allocated<-system2('scontrol',c('show','hostnames',shQuote(nodes)),stdout=TRUE)
 short<-function(x)sub('[.].*$','',x)
 if(!short(unname(Sys.info()[['nodename']]))%in%short(allocated))stop('Current host is not an allocated SLURM compute node')
}
# Offline source install for TRUBA. Never modifies a library outside the explicit new target.
a<-commandArgs(TRUE)
lib<-if(length(a))a[1] else file.path('results','dependencies',paste0('R-',getRversion(),'-',R.version$platform),'library')
reuse<-'--reuse-existing' %in% a
if(as.character(getRversion())!='4.5.2')stop('Run with actual R 4.5.2; Rscript --version alone does not verify the interpreter')
lockfile<-'config/dependencies/source_lock.csv';lock<-read.csv(lockfile,stringsAsFactors=FALSE,check.names=FALSE)
oldlibs<-.libPaths();lib<-normalizePath(lib,mustWork=FALSE)
if(lib%in%normalizePath(oldlibs,mustWork=FALSE))stop('Refusing to modify an existing R library')
dir.create(lib,recursive=TRUE,showWarnings=FALSE)
lockmd5<-unname(tools::md5sum(lockfile));stamp<-file.path(lib,'TSR_DEPENDENCY_LOCK_MD5')
if(file.exists(stamp)&&readLines(stamp,warn=FALSE)[1]!=lockmd5)stop('Target belongs to another lock; choose a new directory')
.libPaths(c(lib,.Library))
sha256<-function(path) {
 if(nzchar(Sys.which('sha256sum')))z<-system2('sha256sum',shQuote(path),stdout=TRUE)
 else if(nzchar(Sys.which('shasum')))z<-system2('shasum',c('-a','256',shQuote(path)),stdout=TRUE)
 else stop('sha256sum or shasum required')
 strsplit(z,'[[:space:]]+')[[1]][1]
}
log<-list()
for(i in seq_len(nrow(lock))) {
 z<-lock[i,];archive<-file.path('config/dependencies',z$source_file)
 if(!file.exists(archive)||sha256(archive)!=z$sha256)stop('Source SHA256 mismatch: ',z$package)
 target<-file.path(lib,z$package);desc<-file.path(target,'DESCRIPTION')
 if(file.exists(desc)) {
  if(read.dcf(desc)[1,'Version']!=z$version)stop('Conflicting target version; refusing overwrite: ',z$package)
  action<-'already_installed'
 } else {
  copied<-FALSE
  if(reuse)for(candidate in oldlibs) {
   from<-file.path(candidate,z$package);dp<-file.path(from,'DESCRIPTION')
   if(!file.exists(dp))next
   dc<-read.dcf(dp);built<-if('Built'%in%colnames(dc))dc[1,'Built'] else ''
   if(dc[1,'Version']==z$version&&grepl('R 4[.]5[.]',built)) {
    if(!file.copy(from,lib,recursive=TRUE))stop('Copy failed: ',z$package)
    copied<-TRUE;break
   }
  }
  if(copied)action<-'copied_exact_version_existing_R4.5_build' else {
   # Optional local macOS validation cache; not used for Linux/TRUBA installation.
   binary<-file.path('config/dependencies/macos_binary',paste0(z$package,'_',z$version,'.tgz'))
   binarylock<-'config/dependencies/macos_binary_lock.csv'
   if(reuse&&Sys.info()[['sysname']]=='Darwin'&&file.exists(binary)&&file.exists(binarylock)) {
    b<-read.csv(binarylock,stringsAsFactors=FALSE);b<-b[b$package==z$package&b$version==z$version,,drop=FALSE]
    if(nrow(b)!=1||sha256(binary)!=b$sha256)stop('Binary SHA256 mismatch')
    install.packages(binary,repos=NULL,type='binary',lib=lib);action<-'official_CRAN_R4.5_arm64_binary'
   } else {
    install.packages(archive,repos=NULL,type='source',lib=lib,Ncpus=1L,INSTALL_opts='--no-multiarch')
    action<-'pinned_CRAN_source'
   }
  }
 }
 # Fresh child validation below prevents a preloaded namespace from hiding a conflict.
 if(!file.exists(desc)||read.dcf(desc)[1,'Version']!=z$version)stop('Installation/version failed: ',z$package)
 log[[i]]<-data.frame(package=z$package,version=z$version,action=action,library=lib)
 write.csv(do.call(rbind,log),file.path(dirname(lib),'installation.csv'),row.names=FALSE)
 cat('INSTALLED',z$package,z$version,action,'\n')
}
writeLines(lockmd5,stamp)
writeLines(capture.output(sessionInfo()),file.path(dirname(lib),'installer_session.txt'))
cat('DEPENDENCIES_INSTALLED: validate in a fresh R 4.5.2 process using this library\n')
