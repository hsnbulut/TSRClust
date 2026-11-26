revision_metrics <- function(truth,pred,main_k) {
  check_labels(truth,pred)
  v1<-legacy$evaluate_partition(truth,pred,main_k)
  v1$MinorityClusterF1_legacy_v1<-v1$MinorityClusterF1
  v1$MinorityDetectionF1_legacy_v1<-v1$MinorityDetectionF1
  v1$MinorityPrecision_legacy_v1<-v1$MinorityPrecision
  mt<-legacy$hierarchical_match(truth,pred,main_k)
  pos<-pred %in% mt$remaining;y<-truth>main_k
  tp<-sum(pos&y);fp<-sum(pos&!y);fn<-sum(!pos&y);tn<-sum(!pos&!y)
  ratio<-legacy$safe_ratio
  v2<-data.frame(MinorityPrecision_all_extra_v2=ratio(tp,tp+fp,if(any(y))0 else NA_real_),
    MinorityRecall_all_extra_v2=ratio(tp,tp+fn,if(any(y))0 else NA_real_),
    MinorityDetectionF1_all_extra_v2=if(any(y))ratio(2*tp,2*tp+fp+fn,0) else NA_real_,
    MacroFalsePositiveReassignment=if(any(truth%in%seq_len(main_k))) mean(pos[truth%in%seq_len(main_k)]) else NA_real_)
  list(metrics=cbind(v1,v2),contingency=table(truth,pred),matching=mt,
       confusion_all_extra_v2=matrix(c(tn,fn,fp,tp),2,dimnames=list(truth=c('other','minority'),prediction=c('other','extra'))),
       confusion_legacy_v1=local({q<-pred %in% mt$matched_minority_pred[mt$matched_minority_pred>0];table(factor(y,levels=c(FALSE,TRUE)),factor(q,levels=c(FALSE,TRUE)))}))
}
