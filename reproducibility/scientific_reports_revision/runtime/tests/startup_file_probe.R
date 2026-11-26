# Capture startup state BEFORE the project bootstrap, then use the existing probe.
a<-commandArgs(TRUE);stopifnot(length(a)==2L)
write.csv(data.frame(R_version=as.character(getRversion()),R_LIBS_USER=Sys.getenv('R_LIBS_USER'),
 TSR_R_LIB=Sys.getenv('TSR_R_LIB'),first_lib=.libPaths()[1L],
 R_ENVIRON=Sys.getenv('R_ENVIRON'),R_ENVIRON_USER=Sys.getenv('R_ENVIRON_USER'),
 R_PROFILE=Sys.getenv('R_PROFILE'),R_PROFILE_USER=Sys.getenv('R_PROFILE_USER')),
 paste0(a[2L],'.paths.csv'),row.names=FALSE)
source('tests/child_library_probe.R')
# raw probe quits the process; normal probe returns its full read-only result.
write.csv(readRDS(a[2L])$checks,paste0(a[2L],'.packages.csv'),row.names=FALSE)
