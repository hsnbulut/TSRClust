fit_revision <- function(row,rep_id) {
  ds<-seed_id(row$dataset_id,rep_id,0);als<-seed_id(row$dataset_id,rep_id,1,row$algorithm_repeat);ns<-seed_id(row$dataset_id,rep_id,2,row$null_repeat)
  seedmeta<-data.frame(replication=rep_id,data_seed=ds,algorithm_seed=als,null_seed=ns,algorithm_version=REVISION_VERSION,metric_version='legacy_v1+all_extra_v2')
  warnings<-character();beg<-proc.time()[3];diagout<-list()
  ans<-tryCatch(withCallingHandlers({
    d<-revision_data(row,ds); p0<-proc.time()[3];x<-if(row$mode=='residual_only')d$x else legacy$robust_standardize(d$x);pre<-proc.time()[3]-p0
    fit0<-proc.time()[3];set.seed(als)
    k<-row$main_k;alpha<-row$requested_alpha;macro<-NULL;micro<-NULL;selection<-NULL
    minpts<-max(5L,row$minpts_multiplier*ncol(x));method<-row$method
    # D2 alpha is prespecified; no labels supplied to its selector.
    if(row$regime=='data_selected')alpha<-.1
    if(method=='TSRClust'&&row$mode=='residual_only') {
      raw<-integer(nrow(x));r0<-seq_len(nrow(x));k<-0L;macrosec<-0
    } else if(method %in% c('TSRClust','TCLUST','Oracle_TCLUST','Trimmed_kmeans')) {
      if(method=='Oracle_TCLUST') {k<-row$main_k+row$q;alpha<-row$requested_alpha}
      if(row$regime=='data_selected') {
        ic<-tclust::tclustIC(x,kk=seq_len(row$kmax),cc=row$restr_fact,alpha=row$requested_alpha,whichIC='MIXMIX',nstart=row$nstart,
          niter1=3,niter2=20,nkeep=5,parallel=FALSE,trace=FALSE,store_x=FALSE)
        v<-ic$MIXMIX;v[!is.finite(v)]<-Inf;if(all(is.infinite(v)))stop('No finite MIXMIX')
        pos<-which(v==min(v),arr.ind=TRUE)[1,];k<-ic$kk[pos[1]];selection<-list(criterion='MIXMIX',values=v,kmax=row$kmax)
      }
      if(method=='Trimmed_kmeans') {
        macro<-legacy$fit_tkmeans_fixed(x,k,alpha,row$nstart,legacy$TSR_CONFIG)
        raw<-legacy$normalise_tkmeans_labels(macro$cluster,k)
      } else {
        macro<-legacy$fit_tclust_fixed(x,k,alpha,row$restr_fact,row$nstart,legacy$TSR_CONFIG);raw<-as.integer(macro$cluster)
      }
      r0<-which(raw==0L);macrosec<-proc.time()[3]-fit0
    }
    if(method=='TSRClust') {
      ref<-if(row$mode=='residual_only')NULL else revision_reference(x,macro,row$reference)
      diagout$reference<-reference_diagnostics(ref);diagout$macro_sizes<-table(raw[raw>0]);diagout$selection<-selection
      set.seed(ns);t0<-proc.time()[3]
      micro<-revision_micro(x[r0,,drop=FALSE],minpts,row$B,row$structure_alpha,row$density,ref,row$novelty,row$rule,row$eps_grid_size,
                            conditional_allowed=row$mode=='residual_only')
      microsec<-proc.time()[3]-t0;pred<-raw;ok<-micro$cluster>0
      pred[r0[ok]]<-k+micro$cluster[ok]
    } else if(method %in% c('TCLUST','Trimmed_kmeans','Oracle_TCLUST')) {pred<-raw;microsec<-0
    } else {
      macrosec<-0;microsec<-0
      fit<-switch(method,
        DBSCAN_full=legacy$fit_full_dbscan(x),HDBSCAN_full=legacy$fit_full_hdbscan(x),
        OTRIMLE_fixedK=legacy$fit_otrimle_fixed(x,G=k),
        tEIGEN_fixedK=legacy$fit_teigen_bic(x,Gs=k),tEIGEN_BIC=legacy$fit_teigen_bic(x,Gs=seq_len(row$kmax)),stop('Unknown method'))
      if(method=='OTRIMLE_fixedK'&&!is.null(fit$code)&&fit$code==0)stop('OTRIMLE non-solution code 0')
      pred<-fit$cluster;raw<-pred;r0<-which(pred==0L);diagout$fit_diagnostics<-fit[setdiff(names(fit),c('fit','cluster'))]
      if(!is.null(fit$selected_G))k<-fit$selected_G
    }
    check_labels(d$truth,pred)
    fitsec<-proc.time()[3]-fit0
    met<-revision_metrics(d$truth,pred,row$main_k)
    diagout$evaluation<-met[setdiff(names(met),'metrics')]
    diagout$truth<-d$truth;diagout$predicted<-pred;diagout$raw_cluster<-raw;diagout$residual_index<-r0;diagout$micro<-micro
    native<-if(method=='TSRClust')length(unique(micro$cluster[micro$cluster>0])) else NA_integer_
    cbind(data.frame(status='OK',error='',estimated_K=k,n_residual=length(r0),n_accepted_micro=native,macro_objective=if(is.null(macro)||is.null(macro$obj))NA_real_ else macro$obj,macro_code=if(is.null(macro)||is.null(macro$code))NA_integer_ else macro$code,
      any_native_micro=if(is.na(native))NA_real_ else as.numeric(native>0),
      residual_assignment_rate=if(method=='TSRClust'&&length(r0))mean(micro$cluster>0) else NA_real_,
      preprocessing_seconds=pre,macro_seconds=macrosec,micro_seconds=microsec,fit_seconds=fitsec,
      end_to_end_seconds=pre+fitsec,selected_eps=if(is.null(micro))NA_real_ else micro$eps,
      global_p=if(is.null(micro))NA_real_ else micro$p_value,rejection_reason=if(is.null(micro))'' else micro$reason),met$metrics)
  },warning=function(w) {warnings<<-c(warnings,conditionMessage(w));invokeRestart('muffleWarning')}),
  error=function(e)data.frame(status='ERROR',error=conditionMessage(e),end_to_end_seconds=proc.time()[3]-beg))
  ans$warnings<-paste(unique(warnings),collapse=' | ');ans$total_attempt_seconds<-proc.time()[3]-beg
  ans$result_object_bytes<-as.numeric(object.size(diagout))
  list(row=cbind(row,seedmeta,ans),diagnostics=diagout)
}
