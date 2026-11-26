# Secondary analysis only. Not included in any Monte Carlo manifest.
source('R/load.R')
a<-commandArgs(TRUE);if(!length(a))stop('Usage: shuttle_audit.R archived | saved_fit.rds | refit')
out<-'results/shuttle';dir.create(out,recursive=TRUE,showWarnings=FALSE)
if(a[1]=='archived') {
 comp<-read.csv(file.path(SOURCE_ROOT,'github_results/real_data_results/shuttle_microcluster_composition.csv'))
 cls<-read.csv(file.path(SOURCE_ROOT,'github_results/real_data_results/class_trimming_rates.csv'))
 cls<-cls[cls$dataset=='Statlog_Shuttle',];b<-cls[cls$local_class=='Bypass',]
 if(!all(comp$dominant_class=='Bypass'&comp$purity==1))stop('Archived composition insufficient for exact counts')
 for(cl in cls$local_class)comp[[paste0('class_',make.names(cl))]]<-if(cl=='Bypass')comp$size else 0
 comp$Bypass_full_percent<-100*comp$size/b$total;comp$Bypass_trimmed_percent<-100*comp$size/b$trimmed
 write.csv(comp,file.path(out,'archived_microcluster_composition.csv'),row.names=FALSE)
 s<-data.frame(source='submitted archived summary (no memberships)',Bypass_total=b$total,Bypass_in_R0=b$trimmed,
  Bypass_recovered=sum(comp$size),n_covering_clusters=nrow(comp),Bypass_full_recall=sum(comp$size)/b$total,
  Bypass_R0_recall=sum(comp$size)/b$trimmed,residual_assigned_fraction=sum(comp$size)/sum(cls$trimmed))
 write.csv(s,file.path(out,'archived_coverage.csv'),row.names=FALSE);print(s)
 quit(save='no')
}
if(a[1]=='refit') {
 # Explicitly invoked refit; run on a compute node for full Shuttle.
 if(Sys.getenv('TSR_SHUTTLE_COMPUTE_CONFIRMED')!='1')stop('Set TSR_SHUTTLE_COMPUTE_CONFIRMED=1 on a compute node')
 if(!requireNamespace('mlbench',quietly=TRUE))stop('mlbench required')
 data('Shuttle',package='mlbench',envir=environment())
 x<-legacy$robust_standardize(as.matrix(Shuttle[,1:9]));truth<-Shuttle[[10]];k<-2L
 set.seed(legacy$TSR_V2$base_seed+100L);macro<-legacy$fit_tclust_fixed(x,k,.05,20,500,legacy$TSR_CONFIG)
 set.seed(legacy$TSR_V2$base_seed+200L);fit<-legacy$fit_tsr_from_macro(x,macro$cluster,k,TRUE,49)
 obj<-list(x=x,truth=truth,raw=macro$cluster,final=fit$cluster,k=k,micro=fit$micro,source='new submitted-seed refit')
 atomic_rds(obj,file.path(out,'refit_memberships.rds'))
} else obj<-readRDS(a[1])
stopifnot(all(c('x','truth','raw','final','k')%in%names(obj)))
x<-as.matrix(obj$x);truth<-as.character(obj$truth);raw<-obj$raw;final<-obj$final;k<-obj$k
check_labels(truth,final);stopifnot(nrow(x)==length(truth),length(raw)==length(truth))
r0<-which(raw==0);ids<-sort(unique(final[final>k]));classes<-sort(unique(truth));bypass<-truth=='Bypass'
ratio<-legacy$safe_ratio
comp<-bind_rows(lapply(ids,function(g) {
 ix<-which(final==g);counts<-vapply(classes,function(cl)sum(truth[ix]==cl),integer(1));z<-data.frame(micro_id=g,size=length(ix),Bypass=sum(bypass[ix]),purity=max(counts)/length(ix),Bypass_full_percent=100*ratio(sum(bypass[ix]),sum(bypass)),Bypass_trimmed_percent=100*ratio(sum(bypass[ix]),sum(bypass[r0])))
 cbind(z,as.data.frame(as.list(setNames(counts,paste0('class_',make.names(classes))))))
}))
write.csv(comp,file.path(out,'microcluster_composition.csv'),row.names=FALSE)
entropy<-function(p){p<-p[p>0];-sum(p*log(p))}
info<-function(y,z) {
 t<-table(y,z);p<-t/sum(t);py<-rowSums(p);pz<-colSums(p);hy<-entropy(py);hz<-entropy(pz)
 hyz<-sum(vapply(seq_len(ncol(p)),function(j)if(pz[j])pz[j]*entropy(p[,j]/pz[j]) else 0,numeric(1)))
 hzy<-sum(vapply(seq_len(nrow(p)),function(i)if(py[i])py[i]*entropy(p[i,]/py[i]) else 0,numeric(1)))
 h<-if(hy==0)1 else max(0,1-hyz/hy);c<-if(hz==0)1 else max(0,1-hzy/hz)
 c(homogeneity=h,completeness=c,V_measure=if(h+c)2*h*c/(h+c) else 0)
}
accepted<-final>k;s<-data.frame(source=if(is.null(obj$source))a[1] else obj$source,n_residual=length(r0),Bypass_total=sum(bypass),Bypass_in_R0=sum(bypass[r0]),Bypass_recovered=sum(bypass&accepted),Bypass_full_recall=ratio(sum(bypass&accepted),sum(bypass)),Bypass_R0_recall=ratio(sum(bypass[r0]&accepted[r0]),sum(bypass[r0])),n_covering_clusters=length(unique(final[bypass&accepted])),residual_assigned_fraction=mean(accepted[r0]),ResidualARI=legacy$adjusted_rand_index(truth[r0],final[r0]))
s<-cbind(s,as.data.frame(as.list(info(truth[r0],final[r0]))));write.csv(s,file.path(out,'coverage_and_residual_metrics.csv'),row.names=FALSE)
# Describe the grid and geometry without selecting/merging by class labels.
ref<-legacy$make_macro_reference(x,raw,k);minpts<-max(5,2*ncol(x));rx<-x[r0,,drop=FALSE]
grid<-legacy$make_eps_grid(rx,minpts);diag<-list();overlap<-list()
for(eps in grid) {
 lab<-dbscan::dbscan(rx,eps=eps,minPts=minpts)$cluster
 z<-legacy$cluster_density_scores(rx,lab,minpts,2,ref,qchisq(.995,ncol(x)));z$eps<-rep(eps,nrow(z));diag[[length(diag)+1]]<-z
 overlap[[length(overlap)+1]]<-data.frame(eps=eps,as.data.frame(table(candidate=lab,true_class=truth[r0])))
}
write.csv(bind_rows(diag),file.path(out,'epsilon_candidate_geometry.csv'),row.names=FALSE)
write.csv(bind_rows(overlap),file.path(out,'epsilon_class_composition_diagnostic.csv'),row.names=FALSE)
if(length(ids)>1) {
 centers<-t(vapply(ids,function(g)colMeans(x[final==g,,drop=FALSE]),numeric(ncol(x))))
 write.csv(as.matrix(dist(centers)),file.path(out,'accepted_centroid_distances.csv'),row.names=TRUE)
}
cat('SHUTTLE_AUDIT_OK; interpret grid connectivity, density and novelty jointly. No label-based merging performed.\n')
