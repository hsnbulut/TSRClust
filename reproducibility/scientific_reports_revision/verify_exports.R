# Verify immutable public scalar exports. No observations or model fits are created.
args<-commandArgs(TRUE);stopifnot(length(args)==1)
D<-normalizePath(args[1]);tmp<-tempfile('tsr_public_exports_');dir.create(tmp);on.exit<-NULL
untar(file.path(D,'results/TSRClust_Tier1_Analysis.tar.gz'),exdir=tmp)
x<-read.csv(file.path(tmp,'replications.csv'));s<-read.csv(file.path(tmp,'summary.csv'));p<-read.csv(file.path(tmp,'paired_differences.csv'))
d<-read.csv(file.path(D,'runtime/config/design.csv'));d<-d[d$tier==1,];t<-read.csv(file.path(D,'runtime/config/tasks_tier1.csv'))
stopifnot(nrow(x)==98800,all(x$status=='OK'),length(unique(x$cell_id))==282,nrow(t)==9880,!anyDuplicated(x[c('cell_id','replication')]))
keys<-unlist(lapply(seq_len(nrow(t)),function(i)paste(t$cell_id[i],seq.int(t$rep_start[i],t$rep_end[i]))))
stopifnot(!anyDuplicated(keys),setequal(keys,paste(x$cell_id,x$replication)))
e<-new.env();sys.source(file.path(D,'runtime/R/common.R'),e)
for(stream in c('data','algorithm','null')){
 repeat_id<-if(stream=='data')rep(1L,nrow(x)) else x[[paste0(stream,'_repeat')]]
 stopifnot(all(x[[paste0(stream,'_seed')]]==mapply(e$seed_id,x$dataset_id,x$replication,repeat_id,MoreArgs=list(stream=match(stream,c('data','algorithm','null'))-1L))))
}
close<-function(a,b)if(is.na(a)||is.na(b))is.na(a)&&is.na(b) else abs(a-b)<1e-10
for(i in seq_len(nrow(s))){
 z<-x[x$cell_id==s$cell_id[i],s$metric[i]];z<-z[is.finite(z)]
 stopifnot(length(z)==s$n_finite[i],close(if(length(z))mean(z) else NA_real_,s$mean[i]))
}
for(i in seq_len(nrow(p))){
 a<-x[x$cell_id==p$cell_id[i],];b<-x[x$cell_id==p$baseline_cell[i],];w<-merge(a,b,by='replication',suffixes=c('_new','_base'))
 z<-w[[paste0(p$metric[i],'_new')]]-w[[paste0(p$metric[i],'_base')]];z<-z[is.finite(z)]
 if(!close(if(length(z))mean(z) else NA_real_,p$mean_difference[i])) {print(p[i,]);print(c(actual=if(length(z))mean(z) else NA_real_,expected=p$mean_difference[i],n=length(z)));stop('Paired export mismatch at row ',i)}
 stopifnot(length(z)==p$n_pairs[i])
}
unlink(tmp,recursive=TRUE)
cat('PUBLIC_EXPORTS_OK 98800 records; 282 configurations; 9880 blocks;',nrow(s),'means;',nrow(p),'paired differences; all three seed streams\n')
