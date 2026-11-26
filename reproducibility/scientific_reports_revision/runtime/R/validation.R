mc_p <- function(score,null_scores) (1+sum(null_scores>=score))/(length(null_scores)+1)
mc_accept <- function(score,null_scores,alpha=.05,rule='legacy_rule') {
  score>0 && mc_p(score,null_scores)<=alpha && (rule!='legacy_rule'||score>max(null_scores))
}
conditional_box <- function(x) {
  n<-nrow(x);p<-ncol(x)
  if(n<2L||any(apply(x,2,function(z) diff(range(z)))==0)) stop('nondegenerate ranges and n>=2 required')
  out<-matrix(0,n,p)
  for(j in seq_len(p)) {
    a<-range(x[,j]); if(diff(a)<=0) stop('degenerate margin')
    out[,j]<-sample(c(a,runif(n-2L,a[1],a[2])))
  }
  out
}
revision_reference <- function(x,macro,kind='legacy_empirical') {
  k<-macro$k
  if(kind=='legacy_empirical') return(legacy$make_macro_reference(x,macro$cluster,k))
  if(kind=='tclust_direct') return(lapply(seq_len(k),function(g) {
    list(center=macro$centers[,g],covariance=macro$cov[,,g])
  }))
  if(kind!='regularized_mcd') stop('unknown reference')
  if(!requireNamespace('robustbase',quietly=TRUE)) stop('robustbase required for regularized_mcd')
  lapply(seq_len(k),function(g) {
    xx<-x[macro$cluster==g,,drop=FALSE]
    if(nrow(xx)<=ncol(xx)+1) stop('MCD reference infeasible: too few macro observations')
    m<-robustbase::covMcd(xx)
    e<-eigen(m$cov,symmetric=TRUE);floor<-max(mean(e$values),1)*1e-4
    list(center=m$center,covariance=e$vectors%*%diag(pmax(e$values,floor),ncol(xx))%*%t(e$vectors))
  })
}
reference_diagnostics <- function(ref) bind_rows(lapply(seq_along(ref),function(g) {
  v<-eigen(ref[[g]]$covariance,symmetric=TRUE,only.values=TRUE)$values
  data.frame(macro=g,min_eigen=min(v),condition_number=max(v)/min(v))
}))
revision_micro <- function(x,min_pts,B=49L,alpha=.05,density=2,ref=NULL,novelty=.995,
                           rule='legacy_rule',grid_size=31L,conditional_allowed=FALSE) {
  if(!rule %in% c('legacy_rule','pvalue_only','recalibrated_rule')) stop('unknown rule')
  if(rule=='recalibrated_rule'&&!conditional_allowed) stop('recalibrated_rule requires independent uniform stage-II null')
  empty<-list(cluster=integer(nrow(x)),eps=NA_real_,p_value=1,score=0,threshold=NA_real_,
              n_candidates=0L,null_scores=numeric(),scan=data.frame(),accepted_sizes=integer())
  if(nrow(x)<min_pts) {empty$reason<-'infeasible_minPts';return(empty)}
  threshold<-qchisq(novelty,ncol(x)); grid<-legacy$make_eps_grid(x,min_pts,grid_size)
  # Save all candidates, including failed geometry filters, for rejection diagnostics.
  scans<-lapply(grid,function(eps) {
    labels<-dbscan::dbscan(x,eps=eps,minPts=min_pts)$cluster
    s<-legacy$cluster_density_scores(x,labels,min_pts,density,ref,threshold)
    s$eps<-rep(eps,nrow(s));s$reason<-ifelse(s$size<min_pts,'size',ifelse(s$density_ratio<density,'density',ifelse(s$novelty_distance<threshold,'novelty','eligible')))
    list(eps=eps,labels=labels,scores=s)
  })
  allscan<-bind_rows(lapply(scans,`[[`,'scores'))
  observed<-lapply(scans,function(z) {z$scores<-z$scores[z$scores$score>0,,drop=FALSE];z})
  observed<-Filter(function(z)nrow(z$scores)>0,observed)
  obs<-legacy$max_candidate_score(observed)
  empty$scan<-allscan
  if(obs<=0||!length(observed)) {empty$reason<-'no_geometric_candidate';return(empty)}
  null<-numeric(B)
  for(b in seq_len(B)) {
    xx<-if(rule=='recalibrated_rule') conditional_box(x) else legacy$simulate_diffuse_noise(x)
    gg<-if(rule=='recalibrated_rule') legacy$make_eps_grid(xx,min_pts,grid_size) else grid
    null[b]<-legacy$max_candidate_score(legacy$scan_dbscan_candidates(xx,min_pts,gg,density,ref,threshold))
  }
  p<-mc_p(obs,null);empty$p_value<-p;empty$score<-obs;empty$threshold<-max(null);empty$null_scores<-null
  empty$n_candidates<-sum(vapply(observed,function(z)nrow(z$scores),integer(1)))
  if(p>alpha) {empty$reason<-'global_p';return(empty)}
  eligible<-lapply(observed,function(z) {
    ok<-vapply(z$scores$score,mc_accept,logical(1),null_scores=null,alpha=alpha,rule=rule)
    z$accepted_clusters<-z$scores$cluster[ok]
    z$quality<-sum(z$scores$size[ok]*z$scores$density_ratio[ok]^2)
    z$coverage<-sum(z$scores$size[ok]); z
  })
  eligible<-Filter(function(z)length(z$accepted_clusters)>0,eligible)
  if(!length(eligible)) {empty$reason<-'candidate_max';return(empty)}
  best<-eligible[[order(vapply(eligible,`[[`,numeric(1),'quality'),vapply(eligible,`[[`,numeric(1),'coverage'),decreasing=TRUE)[1]]]
  labels<-best$labels;labels[!labels %in% best$accepted_clusters]<-0L
  pos<-sort(unique(labels[labels>0]));labels[labels>0]<-match(labels[labels>0],pos)
  empty$cluster<-as.integer(labels);empty$eps<-best$eps;empty$accepted_sizes<-as.integer(table(labels[labels>0]));empty$reason<-'accepted'
  if(nrow(empty$scan)) empty$scan$candidate_p<-vapply(empty$scan$score,mc_p,numeric(1),null_scores=null)
  empty
}
