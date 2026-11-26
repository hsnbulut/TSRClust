source('R/common.R')
dir.create('config',showWarnings=FALSE);dir.create('tables',showWarnings=FALSE)
base<-data.frame(n=500L,p=5L,main_k=2L,separation=6,main_structure='elliptical',main_distribution='normal',main_weights='balanced',residual_total=.1,minority_ratio=.05,q=1L,minority_shape='gaussian',minority_sd=.2,minority_proximity='isolated',noise_type='uniform',mode='end_to_end',nR=50L)
scenario<-function(id,...) {z<-base;for(nm in names(list(...))) z[[nm]]<-list(...)[[nm]];cbind(scenario=id,z)}
ss<-list(); add_s<-function(id,...) {ss[[length(ss)+1L]]<<-scenario(id,...);id}
nulls<-character()
for(noise in c('uniform','correlated_gaussian','heterogeneous','heavy_tail','bridge','correlated_sheet')) for(p in c(2L,5L)) for(r in c(.1,.2)) nulls<-c(nulls,add_s(paste('N',noise,p,r,sep='_'),p=p,residual_total=r,minority_ratio=0,q=0L,noise_type=noise))
cond<-character();for(p in c(2L,5L)) for(nr in c(50L,100L)) cond<-c(cond,add_s(paste('COND',p,nr,sep='_'),p=p,mode='residual_only',nR=nr,minority_ratio=0,q=0L))
power<-character();for(p in c(2L,5L)) for(m in c(.025,.05)) power<-c(power,add_s(paste('POWER',p,m,sep='_'),p=p,minority_ratio=m))
bcases<-character();for(m in c(.025,.05)) for(prox in c('isolated','adjacent','overlapping')) bcases<-c(bcases,add_s(paste('B',m,prox,sep='_'),minority_ratio=m,minority_proximity=prox))
positive<-add_s('positive');negative<-add_s('negative',q=0L,minority_ratio=0)
weak<-add_s('weak',separation=4,minority_sd=.55)
high<-add_s('high_p25',p=25L,n=1000L)
multi<-add_s('two_minorities',q=2L)
tier2null<-character();for(noise in c('uniform','correlated_gaussian','heavy_tail')) tier2null<-c(tier2null,add_s(paste0('T2_',noise),p=10L,noise_type=noise,q=0L,minority_ratio=0))
scenarios<-bind_rows(ss); settings<-setdiff(names(scenarios),'scenario')
keys<-apply(scenarios[settings],1,paste,collapse='|');scenarios$dataset_id<-match(keys,unique(keys))
write.csv(scenarios,'config/scenarios.csv',row.names=FALSE)
basecfg<-data.frame(tier=1L,method='TSRClust',rule='legacy_rule',reference='legacy_empirical',novelty=.995,density=2,minpts_multiplier=2L,B=49L,nstart=500L,regime='fixed_macro',kmax=5L,algorithm_repeat=1L,null_repeat=1L,oracle=FALSE,structure_alpha=.05,restr_fact=20,eps_grid_size=31L)
cells<-list(); add<-function(exp,sc,reps,...) {z<-basecfg;for(nm in names(list(...)))z[[nm]]<-list(...)[[nm]];cells[[length(cells)+1L]]<<-cbind(experiment=exp,scenario=sc,reps=reps,z)}
for(sc in nulls)for(rule in c('legacy_rule','pvalue_only'))add('A-END',sc,1000L,rule=rule)
for(sc in cond)for(rule in c('legacy_rule','pvalue_only','recalibrated_rule'))add('A-COND',sc,1000L,rule=rule)
for(sc in power)for(rule in c('legacy_rule','pvalue_only'))add('A-POWER',sc,200L,rule=rule)
for(sc in bcases)for(nov in c(.95,.975,.99,.995,.999))for(ref in c('legacy_empirical','tclust_direct','regularized_mcd'))add('B',sc,200L,novelty=nov,reference=ref)
for(sc in c(positive,negative,weak,high)) {
 add('C',sc,200L)
 for(v in c(1.5,2.5,3))add('C',sc,200L,density=v)
 for(v in c(1L,3L,4L))add('C',sc,200L,minpts_multiplier=v)
 for(v in c(99L,199L))add('C',sc,200L,B=v)
 for(v in c(.95,.975,.99,.999))add('C',sc,200L,novelty=v)
 for(de in c(1.5,3))for(mp in c(1L,4L))add('C',sc,200L,density=de,minpts_multiplier=mp)
}
for(sc in c(negative,positive,multi,weak)) {
 for(m in c('TSRClust','TCLUST','Trimmed_kmeans','DBSCAN_full','HDBSCAN_full','OTRIMLE_fixedK','tEIGEN_fixedK','Oracle_TCLUST'))add('D1-E',sc,100L,method=m,regime='controlled_K',oracle=m=='Oracle_TCLUST')
 for(kmax in c(3L,5L,7L))for(m in c('TSRClust','tEIGEN_BIC'))add('D2-E',sc,100L,method=m,regime='data_selected',kmax=kmax)
}
for(sc in tier2null)for(B in c(49L,99L,199L))for(rule in c('legacy_rule','pvalue_only'))add('A-END',sc,1000L,tier=2L,B=B,rule=rule)
for(sc in c(positive,negative,weak))for(ar in 1:2)for(nr in 1:2)add('C-RNG',sc,100L,tier=2L,algorithm_repeat=ar,null_repeat=nr)
cfg<-bind_rows(cells);cfg$cell_id<-seq_len(nrow(cfg))
design<-merge(cfg,scenarios,by='scenario',sort=FALSE);design<-design[order(design$cell_id),];rownames(design)<-NULL
design$requested_alpha<-ifelse(design$regime=='data_selected',.1,ifelse(design$oracle,design$residual_total-design$minority_ratio,design$residual_total))
design$K_information<-ifelse(design$oracle,'oracle_total',ifelse(design$regime=='data_selected','unlabeled_IC',ifelse(design$method%in%c('DBSCAN_full','HDBSCAN_full'),'none','true_dominant_K')))
# minPts changes in p=2 may intentionally coincide because max(5,p) is floor-limited.
write.csv(design,'config/design.csv',row.names=FALSE)
make_tasks<-function(d,pilot=FALSE) {
 xs<-list()
 for(i in seq_len(nrow(d)))for(start in seq(1L,if(pilot)2L else d$reps[i],by=if(pilot)2L else 10L))xs[[length(xs)+1L]]<-data.frame(cell_id=d$cell_id[i],rep_start=start,rep_end=min(start+if(pilot)1L else 9L,if(pilot)2L else d$reps[i]))
 out<-bind_rows(xs);out$task_id<-seq_len(nrow(out));out
}
# Representative and bounded pilot, including unavailable optional packages as explicit failures.
pick<-c(which(design$experiment=='A-END'&design$p==2&design$noise_type=='uniform')[1:2],
 which(design$experiment=='A-COND')[1:3],which(design$experiment=='B'&design$novelty==.995)[1:3],
 which(design$experiment=='C'&design$scenario %in% c('positive','high_p25')&design$density==2&design$minpts_multiplier==2&design$B==49&design$novelty==.995),
 which(design$experiment=='D1-E'&design$scenario=='positive'),
 which(design$experiment=='D2-E'&design$scenario=='positive'&design$kmax==3))
pick<-unique(pick[!is.na(pick)])
write.csv(make_tasks(design[pick,],TRUE),'config/tasks_pilot.csv',row.names=FALSE)
for(tier in 1:2)write.csv(make_tasks(design[design$tier==tier,]),paste0('config/tasks_tier',tier,'.csv'),row.names=FALSE)
design$outer_fits<-design$reps
sel<-design$method=='TSRClust'&design$regime=='data_selected'
design$macro_fit_upper<-ifelse(design$mode=='residual_only'|design$method%in%c('DBSCAN_full','HDBSCAN_full','tEIGEN_fixedK','tEIGEN_BIC','OTRIMLE_fixedK'),0,ifelse(sel,design$kmax+1,1))*design$reps
# teigen: candidate G count only; covariance models and iterative fitting add work.
design$teigen_G_candidates<-ifelse(design$method=='tEIGEN_BIC',design$kmax,ifelse(design$method=='tEIGEN_fixedK',1,0))*design$reps
 design$null_upper<-ifelse(design$method=='TSRClust',design$B,0)*design$reps
work<-aggregate(cbind(outer_fits,macro_fit_upper,teigen_G_candidates,null_upper)~tier+experiment,design,sum)
work$cells<-as.integer(table(paste(design$tier,design$experiment))[paste(work$tier,work$experiment)])
write.csv(work,'tables/workload.csv',row.names=FALSE)
write.csv(design[c('cell_id','scenario','experiment','tier','reps','outer_fits','macro_fit_upper','teigen_G_candidates','null_upper')],'tables/workload_by_cell.csv',row.names=FALSE)
cat('DESIGN_OK:',nrow(design),'configuration cells;',nrow(make_tasks(design[pick,],TRUE)),'pilot tasks\n');print(work)
