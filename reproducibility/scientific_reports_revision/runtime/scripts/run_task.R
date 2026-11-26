# Fail before any native package can be loaded by an unintended R launcher.
if(as.character(getRversion())!='4.5.2')stop('Worker requires actual R 4.5.2; check the interpreter and R module')
source('R/library_environment.R')
initialize_revision_library(required=TRUE)
dep_lock<-read.csv('config/dependencies/source_lock.csv',stringsAsFactors=FALSE)
for(i in seq_len(nrow(dep_lock))) {
 p<-dep_lock$package[i];desc<-tryCatch(utils::packageDescription(p),error=function(e)NULL)
 if(is.null(desc)||is.na(desc$Version)||desc$Version!=dep_lock$version[i])stop('Pinned dependency unavailable/version mismatch: ',p)
 if(identical(desc$NeedsCompilation,'yes')&&!grepl('R 4[.]5[.]',desc$Built))stop('Native dependency is not an R 4.5 build: ',p)
}
source('R/load.R')
a<-commandArgs(TRUE);if(length(a)<2)stop('Usage: run_task.R config/tasks_pilot.csv task_id [output_dir]')
manifest<-a[1];taskid<-as.integer(a[2]);outdir<-if(length(a)>2)a[3] else 'results/pilot'
if(!grepl('pilot',basename(manifest))&&Sys.getenv('TSR_FULL_APPROVED')!='REVIEWED_AND_APPROVED')stop('Full execution requires implementation/pilot review; set TSR_FULL_APPROVED only after approval')
tasks<-read.csv(manifest);task<-tasks[tasks$task_id==taskid,,drop=FALSE];if(nrow(task)!=1)stop('Invalid array task index')
design<-read.csv('config/design.csv',stringsAsFactors=FALSE);row<-design[design$cell_id==task$cell_id,,drop=FALSE];if(nrow(row)!=1)stop('Invalid cell')
# Fail closed if the frozen implementation or submitted source has drifted.
pins<-read.csv('config/runtime_md5.csv',stringsAsFactors=FALSE)
actual<-unname(tools::md5sum(pins$path));if(anyNA(actual)||!identical(actual,pins$md5))stop('Frozen source hash mismatch; review and regenerate runtime_md5.csv')
identity<-list(manifest_md5=unname(tools::md5sum(manifest)),design_md5=unname(tools::md5sum('config/design.csv')),
 runtime_md5=unname(tools::md5sum('config/runtime_md5.csv')),version=REVISION_VERSION)
path<-file.path(outdir,sprintf('task_%06d.rds',taskid));dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
lock<-paste0(path,'.lock');if(!dir.create(lock,showWarnings=FALSE))stop('Task locked; verify no running writer before removing ',lock)
on.exit_cleanup<-NULL
tryCatch({
 previous<-if(file.exists(path))readRDS(path) else NULL
 if(!is.null(previous)&&!identical(previous$identity,identity))stop('Output belongs to another frozen design/source; use a new output directory')
 results<-if(is.null(previous))list() else previous$results
 for(rep_id in seq.int(task$rep_start,task$rep_end)) {
  key<-as.character(rep_id)
  if(!is.null(results[[key]])&&results[[key]]$row$status=='OK')next
  attempts<-if(is.null(results[[key]]))list() else c(results[[key]]$previous_attempts,list(results[[key]][c('row','diagnostics')]))
  cpu_start<-proc.time()
  z<-fit_revision(row,rep_id)
  delta<-proc.time()-cpu_start
  z$row$replication_cpu_seconds<-unname(delta[['user.self']]+delta[['sys.self']])
  z$row$replication_elapsed_seconds<-unname(delta[['elapsed']])
  z$row$slurm_job_id<-Sys.getenv('SLURM_JOB_ID')
  z$row$slurm_array_task_id<-Sys.getenv('SLURM_ARRAY_TASK_ID')
  z$previous_attempts<-attempts;results[[key]]<-z
  atomic_rds(list(identity=identity,task=task,results=results,session=capture.output(sessionInfo()),rng=RNGkind()),path)
 }
 cat('TASK_FINISHED',taskid,'rows',length(results),'failed',sum(vapply(results,function(z)z$row$status!='OK',logical(1))),'\n')
},finally=unlink(lock,recursive=TRUE))
