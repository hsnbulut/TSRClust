source('R/load.R')
rscript<-Sys.getenv('TSR_RSCRIPT',normalizePath('scripts/tsr_rscript.sh'))
rscript_args<-character()
probe<-file.path(tempdir(),'tsr_workflow_checks');dir.create(probe,showWarnings=FALSE)
run<-function(script,args) {
 z<-system2(rscript,c(rscript_args,shQuote(script),vapply(args,shQuote,character(1))),stdout=TRUE,stderr=TRUE)
 code<-attr(z,'status');if(!is.null(code)&&code!=0)stop(paste(z,collapse='\n'));z
}
run('scripts/run_task.R',c('config/tasks_pilot.csv','1',probe))
f<-file.path(probe,'task_000001.rds');hash<-tools::md5sum(f);mtime<-file.info(f)$mtime
run('scripts/run_task.R',c('config/tasks_pilot.csv','1',probe))
stopifnot(identical(hash,tools::md5sum(f)),identical(mtime,file.info(f)$mtime))
run('scripts/check_complete.R',c('config/tasks_pilot.csv',probe));z<-read.csv(file.path(probe,'completeness.csv'))
stopifnot(sum(z$successful)==2,sum(z$missing)==38,sum(z$failed)==0)
# Controlled checkpoint failure injection in a temporary test result, never in scientific output.
obj<-readRDS(f);obj$results[['2']]$row$status<-'ERROR';obj$results[['2']]$row$error<-'controlled recovery test'
saveRDS(obj,f)
run('scripts/run_task.R',c('config/tasks_pilot.csv','1',probe))
recovered<-readRDS(f)
stopifnot(recovered$results[['2']]$row$status=='OK',length(recovered$results[['2']]$previous_attempts)==1,
 recovered$results[['2']]$previous_attempts[[1]]$row$error=='controlled recovery test')
# Full run protection is exercised without fitting.
Sys.unsetenv('TSR_FULL_APPROVED')
z<-suppressWarnings(system2(rscript,c(rscript_args,'scripts/run_task.R','config/tasks_tier1.csv','1',shQuote(probe)),stdout=TRUE,stderr=TRUE))
stopifnot(attr(z,'status')!=0,any(grepl('Full execution requires',z)))
write.csv(data.frame(check=c('atomic checkpoint/resume preserves successful output','missing tasks exposed','failed checkpoint recovers and preserves attempt','full execution gated'),status='PASS'),file.path(Sys.getenv('TSR_ACCEPTANCE_RUN_DIR','reports'),'workflow_checks.csv'),row.names=FALSE)
unlink(probe,recursive=TRUE);cat('WORKFLOW_CHECKS_OK\n')
