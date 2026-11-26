# Evaluate only the Linux allocation guard with controlled environment functions.
# This never invokes package installation, an actual scheduler, or remote commands.
guard<-parse('scripts/install_dependencies.R')[[1L]]
probe<-function(host,job,nodes) {
 e<-new.env(parent=baseenv())
 e$Sys.info<-function()c(sysname='Linux',nodename=host)
 e$Sys.getenv<-function(name,unset='')switch(name,SLURM_JOB_ID=job,SLURM_JOB_NODELIST=nodes,SLURM_NODELIST=nodes,unset)
 e$Sys.which<-function(...)'/stub/scontrol'
 e$system2<-function(...)'orfoz001'
 tryCatch({eval(guard,e);'PASS_COMPUTE'},error=function(e)conditionMessage(e))
}
stopifnot(grepl('requires a SLURM compute-node allocation',probe('arf-ui1','','')))
stopifnot(grepl('not an allocated SLURM compute node',probe('arf-ui1','123','orfoz001')))
stopifnot(probe('orfoz001.yonetim','123','orfoz001')=='PASS_COMPUTE')
write.csv(data.frame(check=c('Linux login without allocation is rejected','login host with salloc ID is rejected','allocated compute hostname passes'),status='PASS'),'reports/truba_compute_guard_checks.csv',row.names=FALSE)
cat('COMPUTE_GUARD_CONTRACTS_OK 3\n')
