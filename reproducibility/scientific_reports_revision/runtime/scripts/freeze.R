source('config/paths.R')
files<-c(list.files('R',full.names=TRUE,pattern='[.]R$'),list.files('scripts',full.names=TRUE,pattern='[.](R|py|sh)$'),list.files('tests',full.names=TRUE,pattern='[.](R|py)$'),list.files('slurm',full.names=TRUE,pattern='[.](sh|sbatch)$'),list.files('config',full.names=TRUE,pattern='tasks_.*[.]csv$'),
 'config/paths.R','config/design.csv','config/scientific_integrity.csv','config/truba.env','config/truba.env.example','config/dependencies/source_lock.csv','config/dependencies/macos_binary_lock.csv',list.files(file.path(SOURCE_ROOT,'reproducibility/core'),full.names=TRUE,pattern='[.]R$'),
 list.files(file.path(SOURCE_ROOT,'reproducibility/competitors'),full.names=TRUE,pattern='[.]R$'),
 list.files(file.path(SOURCE_ROOT,'reproducibility/evaluation'),full.names=TRUE,pattern='[.]R$'),
 list.files(file.path(SOURCE_ROOT,'reproducibility/config'),full.names=TRUE,pattern='[.]R$'))
write.csv(data.frame(path=files,md5=unname(tools::md5sum(files))),'config/runtime_md5.csv',row.names=FALSE)
cat('FROZEN',length(files),'files\n')
