# Independent read-only final assessment. No fits, generated observations or installs.
stopifnot(as.character(getRversion())=='4.5.2')
arg<-grep('^--file=',commandArgs(),value=TRUE)[1];F<-dirname(normalizePath(sub('^--file=','',arg)))
C3<-dirname(F);R<-dirname(dirname(C3));P<-dirname(R)
archive<-Sys.getenv('TSR_TIER1_ARCHIVE')
sha<-function(p)strsplit(system2('shasum',c('-a','256',shQuote(p)),stdout=TRUE),'[[:space:]]+')[[1]][1]
inv<-read.csv(file.path(F,'evidence/scientific_inputs.csv'))
for(i in seq_len(nrow(inv))){p<-if(nzchar(inv$project_relative[i]))file.path(P,inv$project_relative[i]) else archive;stopifnot(sha(p)==inv$sha256[i])}
wr<-function(x,n,where='tables')write.csv(x,file.path(F,where,n),row.names=FALSE,na='NA')
bind<-function(x)do.call(rbind,x)
numeric_checks<-list()
integer_count<-function(v,name,cap){
 stopifnot(all(is.finite(v)),all(v>=-1e-7),all(v<=cap+1e-7));res<-max(abs(v-round(v)));stopifnot(res<1e-7)
 numeric_checks[[length(numeric_checks)+1]]<<-data.frame(quantity=name,n_values=length(v),max_integer_residual=res,tolerance=1e-7,status='PASS')
 as.integer(round(v))
}
# Independent origin-flow enumeration: NEVER calls the previous origin_bounds.R.
# Enumerate macro E counts; enforce nonnegative M/U/N flow and minimum size of W.
# This is a relaxation of full matching/partition information, so bounds are valid
# outer identification bounds, not claimed sharp over all original fitted labels.
flow_bounds<-function(M,U,N,tp,fv,mp,fe,nz,z,s,vcount,directE=NA_integer_){
 A<-M-mp;other<-fe-fv
 stopifnot(other>=0,s>=0,vcount>=0,tp>=vcount,if(s==0)other==0 else TRUE)
 xe<-seq.int(0L,min(A,fe))
 mz<-A-xe;uz<-z-nz-mz;np<-N-nz-(fe-xe)
 minUW<-max(0L,s-other)
 ok<-mz>=0 & uz>=0 & np>=0 & uz<=U-tp-minUW
 if(is.finite(directE))ok<-ok&xe==directE
 xe<-xe[ok];uz<-uz[ok];np<-np[ok]
 vl<-pmax(0L,xe-other);vu<-pmin(fv,xe)
 ok<-vl<=vu
 xe<-xe[ok];uz<-uz[ok];np<-np[ok];vl<-vl[ok];vu<-vu[ok]
 stopifnot(length(xe)>0)
 lv<-min(vl);uv<-max(vu);le<-min(xe);ue<-max(xe)
 # Construct an aggregate nonnegative integer flow witness for both V endpoints.
 # Witnesses prove compatibility with the stated relaxation,not actual labels.
 for(xv in unique(c(lv,uv))){
  j<-which(vl<=xv & vu>=xv)[1];e<-xe[j];u0<-uz[j];uw<-minUW
  tab<-rbind(M=c(mp,xv,e-xv,A-e),U=c(U-tp-u0-uw,tp,uw,u0),N=c(np[j],fv-xv,other-e+xv,nz))
  stopifnot(all(tab>=0),all(tab==round(tab)),all(rowSums(tab)==c(M,U,N)),sum(tab[,4])==z,sum(tab[c(1,3),2])==fv,sum(tab[c(1,3),2:3])==fe,sum(tab[,3])>=s)
 }
 c(LV=lv,UV=uv,LE=le,UE=ue,n_feasible_macro_E=length(xe))
}
expected_old<-readRDS(file.path(C3,'evidence/identified_aggregate_counts.rds'))
out<-list();cohort_summary<-list();thresholds<-list();collapse<-list();cohort_rates<-list();exports<-list();codes<-list();export_text<-list();headline<-list();witness_n<-0L
audit_block<-function(a,family,cohort,source,direct=FALSE){
 stopifnot(all(a$status=='OK'),!anyDuplicated(a[c('cell_id','replication')]))
 n<-a$n;M<-n-round(n*a$residual_total);U<-round(n*a$minority_ratio);N<-n-M-U
 stopifnot(all(M>0),all(U>=0),all(N>=0),all(U[a$q==0]==0),all(a$main_k>=2))
 vcount<-a$QHat-a$SpuriousExtraClusters;stopifnot(all(vcount>=0),all(vcount<=a$q))
 tp<-integer_count(ifelse(U>0,a$MinorityRecall*U,0),paste(cohort,'TP from recall'),n)
 has<-vcount>0;stopifnot(all((tp>0)==has),all(a$MinorityPrecision[has]>0),all(is.na(a$MinorityClusterF1[U==0])))
 vsize<-rep(0,nrow(a));vsize[has]<-tp[has]/a$MinorityPrecision[has];vsize<-integer_count(vsize,paste(cohort,'V union size from precision'),n)
 fv<-vsize-tp;fn<-U-tp;stopifnot(all(fv>=0),all(fn>=0))
 mp<-integer_count(M*a$MajorRetention,paste(cohort,'macro P'),M)
 fe<-integer_count((M+N)*a$FalseExtraAssignmentRate,paste(cohort,'nonminority E'),M+N)
 nz<-integer_count(ifelse(N>0,N*a$NoiseRecall,0),paste(cohort,'noise zero'),N)
 z<-integer_count(if(direct)a$n_residual else n*a$estimated_noise_proportion,paste(cohort,'hard-label zero total'),n)
 de<-if(direct)integer_count(M*a$MacroFalsePositiveReassignment,paste(cohort,'direct macro E'),M) else rep(NA_integer_,nrow(a))
 yes<-U>0;stopifnot(all(abs(a$MinorityDetectionF1[yes]-2*tp[yes]/(2*tp[yes]+fv[yes]+fn[yes]))<1e-11))
 b<-t(vapply(seq_len(nrow(a)),function(i)flow_bounds(M[i],U[i],N[i],tp[i],fv[i],mp[i],fe[i],nz[i],z[i],a$SpuriousExtraClusters[i],vcount[i],de[i]),numeric(5)))
 # Independently derive and compare the old analytic bound expressions.
 A<-M-mp;wmin<-pmax(0,a$SpuriousExtraClusters-(fe-fv))
 le<-pmax(0,fe-(N-nz),A+nz-z);ue<-pmin(A,fe,A+nz-z+fn-wmin)
 if(direct)le<-ue<-de
 lv<-pmax(0,fv-(N-nz),le-(fe-fv));uv<-pmin(fv,ue)
 le<-pmax(le,lv);ue<-pmin(ue,uv+fe-fv);lv<-pmax(lv,le-(fe-fv));uv<-pmin(uv,ue)
 stopifnot(all(b[,'LV']==lv),all(b[,'UV']==uv),all(b[,'LE']==le),all(b[,'UE']==ue))
 old<-expected_old[expected_old$cohort==cohort,];idx<-match(paste(a$cell_id,a$replication),paste(old$cell_id,old$replication));stopifnot(!anyNA(idx),nrow(old)==nrow(a))
 for(pair in list(c('TP_V','tp'),c('FP_V','fv'),c('FN_V','fn'),c('macro_P','mp'),c('FP_E','fe'),c('noise_0','nz'),c('predicted_0','z'))){stopifnot(all(old[[pair[1]]][idx]==get(pair[2])))}
 for(pair in list(c('macro_V_lower_count','LV'),c('macro_V_upper_count','UV'),c('macro_E_lower_count','LE'),c('macro_E_upper_count','UE')))stopifnot(all(old[[pair[1]]][idx]==b[,pair[2]]))
 stopifnot(all(b[,'LV']>=0),all(b[,'UV']<=b[,'UE']),all(b[,'UE']<=M-mp))
 d<-data.frame(source_family=family,cohort=cohort,cell_id=a$cell_id,replication=a$replication,n=n,n_macro=M,n_minority=U,n_noise=N,TP_V=tp,FP_V=fv,FN_V=fn,V_union_size=vsize,macro_P=mp,FP_E=fe,noise_zero=nz,predicted_zero=z,V_labels=vcount,W_labels=a$SpuriousExtraClusters,macro_V_lower_count=b[,'LV'],macro_V_upper_count=b[,'UV'],macro_E_lower_count=b[,'LE'],macro_E_upper_count=b[,'UE'],lower_rate=b[,'LV']/M,upper_rate=b[,'UV']/M,lower_pct=100*b[,'LV']/M,upper_pct=100*b[,'UV']/M,macro_share_V_upper_pct=ifelse(vsize>0,100*b[,'UV']/vsize,NA_real_),MajorARI=a$MajorARI,MinorityClusterF1=a$MinorityClusterF1,MajorRetention=a$MajorRetention,source=source)
 out[[length(out)+1]]<<-d
 for(target in c('V','E')){
  l<-b[,paste0('L',target)];u<-b[,paste0('U',target)]
  cohort_summary[[length(cohort_summary)+1]]<<-data.frame(source_family=family,cohort=cohort,target=paste0('MacroTo',target),n_fits=nrow(a),n_rate_finite=nrow(a),n_identified=sum(l==u),n_nonunique=sum(l!=u),n_exact_zero=sum(u==0),n_definitely_positive=sum(l>0),n_possibly_positive=sum(u>0),sum_lower_count=sum(l),sum_upper_count=sum(u),mean_lower_rate=mean(l/M),mean_upper_rate=mean(u/M),mean_lower_pct=100*mean(l/M),mean_upper_pct=100*mean(u/M),total_macro_exposures=sum(M),pooled_lower_proportion=sum(l)/sum(M),pooled_upper_proportion=sum(u)/sum(M),max_individual_lower_pct=max(100*l/M),max_individual_upper_pct=max(100*u/M),denominator='All true macro observations per fit;cohort mean averages fitted-record fractions',units='Rates as fractions and explicitly labelled percent;counts summed over repeated simulation observations',interpretation='Identification bounds conditional on saved records;not confidence intervals or population guarantees;valid outer constraints,not claimed globally sharp',source=source)
  for(t in c(0,.001,.005,.01,.05,.10))thresholds[[length(thresholds)+1]]<<-data.frame(source_family=family,cohort=cohort,target=paste0('MacroTo',target),threshold_proportion=t,threshold_percent=100*t,comparison='strictly greater than threshold',n_fits=nrow(a),n_definitely_exceed=sum(l/M>t),n_possibly_exceed=sum(u/M>t),lower_frequency=sum(l/M>t)/nrow(a),outer_upper_frequency=sum(u/M>t)/nrow(a),interpretation='Worst-case outer count of fitted records;not statistical significance or confidence bounds;individual attainability need not imply joint full-matching attainability',source=source)
 }
 for(metric in c('MajorARI','MinorityClusterF1','MajorRetention','MinorityPrecision','MinorityRecall','MinorityDetectionF1','NoiseRecall')){
  val<-a[[metric]];ok<-is.finite(val);cohort_rates[[length(cohort_rates)+1]]<<-data.frame(source_family=family,cohort=cohort,metric=metric,n_fits=nrow(a),n_finite=sum(ok),mean=if(any(ok))mean(val[ok]) else NA_real_,source=source)
 }
 G<-if(direct)a$estimated_K else a$selected_g;stopifnot(all(G==a$main_k));case<-G==2&a$QHat==1&mp>0
 # Exact deterministic original macro group cardinalities,without generating X.
 largest<-vapply(seq_len(nrow(a)),function(i){k<-a$main_k[i];w<-if(a$main_weights[i]=='balanced')rep(1/k,k) else if(k==2)c(.75,.25) else if(k==4)c(.55,.25,.15,.05) else rev(seq_len(k))/sum(seq_len(k));w<-w/sum(w);v<-M[i]*w;cnt<-floor(v);rem<-M[i]-sum(cnt);if(rem>0)cnt[order(v-cnt,decreasing=TRUE)[seq_len(rem)]]<-cnt[order(v-cnt,decreasing=TRUE)[seq_len(rem)]]+1;max(cnt)},numeric(1))
 collapse[[length(collapse)+1]]<<-data.frame(source_family=family,cohort=cohort,n_fits=nrow(a),n_single_macro_matched_component=sum(case),fraction_single_component=mean(case),n_merging_both_macro_groups_certified=sum(case & mp>largest),mean_macro_union_retention_in_subset=if(any(case))mean(a$MajorRetention[case]) else NA,mean_MajorARI_in_subset=if(any(case))mean(a$MajorARI[case]) else NA,mean_MinorityClusterF1_in_subset=if(any(case))mean(a$MinorityClusterF1[case]) else NA,proof='G2,Ecount1,positive macroP imply exactly one macro-matched label;macroP>largest true macro group certifies members from both true groups',source=source)
 witness_n<<-witness_n+nrow(a)
}
work<-tempfile('c3_final_');dir.create(work)
main<-function(){
 on.exit(unlink(work,recursive=TRUE))
 for(st in c('architecture','stress','negative_control')){
  path<-file.path(P,'Scientific Reports Manuscript/analysis/v2/results/final_sim_combined',paste0(st,'_replications.rds'));all<-readRDS(path);a<-all[all$method=='OTRIMLE_fixedK',]
  gzip<-file.path(P,'Scientific Reports Manuscript/github_results/raw_results',paste0(st,'_replications.csv.gz'));ex<-read.csv(gzfile(gzip));ex<-ex[ex$method=='OTRIMLE_fixedK',]
  stopifnot(identical(names(a),names(ex)),nrow(a)==nrow(ex))
  for(nm in names(a)){
   av<-a[[nm]];bv<-ex[[nm]]
   if(is.character(av)){
    if(!identical(is.na(av),is.na(bv)))export_text[[length(export_text)+1]]<<-data.frame(study=st,field=nm,RDS_blank_count=sum(av=='',na.rm=TRUE),RDS_missing_count=sum(is.na(av)),CSV_missing_count=sum(is.na(bv)),normalization='Only for comparison:blank character strings and serialized missing text treated alike;no input files changed')
    av[av=='']<-NA_character_;bv<-as.character(bv);bv[bv=='']<-NA_character_
   }
   stopifnot(identical(is.na(av),is.na(bv)))
   if(is.numeric(av))stopifnot(all(abs(av-bv)<1e-9,na.rm=TRUE)) else stopifnot(all(as.character(av)==as.character(bv),na.rm=TRUE))
  }
  stopifnot(all(a$status=='OK'),all(a$otrimle_code%in%c(1,2)))
  exports[[length(exports)+1]]<<-data.frame(study=st,n_otrimle_rows=nrow(a),n_OK=sum(a$status=='OK'),n_code1=sum(a$otrimle_code==1),n_code2=sum(a$otrimle_code==2),n_code0=sum(a$otrimle_code==0),n_selected_G2=sum(a$selected_g==2),n_selected_G_other=sum(a$selected_g!=2),all_export_fields_checked=ncol(a),RDS_and_gzip_agree_after_blank_text_normalization=TRUE)
  groups<-if(st=='architecture')list('Original architecture M1-M4'=a[a$q>0,],'Original architecture N0'=a[a$q==0,]) else setNames(list(a),paste('Original',st))
  for(co in names(groups)){
   z<-groups[[co]];codes[[length(codes)+1]]<<-data.frame(cohort=co,n_fits=nrow(z),n_code1=sum(z$otrimle_code==1),n_code2=sum(z$otrimle_code==2),inclusion='All originally recorded OK fits retained,including iteration-limit code1;not convergence-only')
   audit_block(z,'Original submitted',co,path)
  }
 }
 untar(archive,files='replications.csv',exdir=work);tier<-read.csv(file.path(work,'replications.csv'));stopifnot(nrow(tier)==98800,all(tier$status=='OK'));a<-tier[tier$method=='OTRIMLE_fixedK',];stopifnot(nrow(a)==400,all(a$experiment=='D1-E'),all(a$MacroFalsePositiveReassignment==0))
 for(sc in unique(a$scenario))audit_block(a[a$scenario==sc,],'Tier1 revision',paste('Tier1',sc),'Tier1 archive replications.csv;D1-E;distinct from originals',TRUE)
 all<-bind(out);saveRDS(all,file.path(F,'evidence/independent_fit_bounds.rds'),compress='xz')
 wr(bind(cohort_summary),'final_macro_promotion_bounds.csv',where='.')
 wr(bind(thresholds),'individual_fit_threshold_diagnostics.csv');wr(bind(cohort_rates),'cohort_metric_means.csv');wr(bind(collapse),'macro_component_collapse.csv');wr(bind(exports),'original_export_checks.csv','evidence');wr(bind(codes),'convergence_inclusion.csv','evidence');wr(bind(numeric_checks),'integer_checks.csv','evidence');wr(bind(export_text),'export_text_normalization.csv','evidence')
 wr(all[all$lower_rate!=all$upper_rate|all$macro_E_lower_count!=all$macro_E_upper_count,],'nonunique_fit_bounds.csv')
 worst<-bind(lapply(unique(all$cohort),function(co){a<-all[all$cohort==co,];a[a$upper_rate==max(a$upper_rate),][1,]}));wr(worst,'maximum_individual_bound_examples.csv')
 cond<-list();worstV<-list()
 for(co in unique(all$cohort)){
  a<-all[all$cohort==co,];ok<-is.finite(a$macro_share_V_upper_pct)
  cond[[length(cond)+1]]<-data.frame(cohort=co,n_fits=nrow(a),n_nonempty_V=sum(ok),max_macro_share_V_upper_pct=if(any(ok))max(a$macro_share_V_upper_pct[ok]) else NA_real_,max_macro_V_upper_count=max(a$macro_V_upper_count),TP_V_sum=sum(a$TP_V),FP_V_sum=sum(a$FP_V),FN_V_sum=sum(a$FN_V),noise_FP_V_lower_count=sum(a$FP_V-a$macro_V_upper_count),noise_FP_V_upper_count=sum(a$FP_V-a$macro_V_lower_count),denominator='macro_share_V uses predicted matched union size,NOT true macro population;cohort totals count repeated exposures',interpretation='Small rMV does not assert negligible relative contamination in every small matched cluster;unavailable for empty V',source=a$source[1])
  if(any(ok)){b<-a[ok,];worstV[[length(worstV)+1]]<-b[b$macro_share_V_upper_pct==max(b$macro_share_V_upper_pct),][1,]}
 }
 wr(bind(cond),'conditional_contamination_diagnostics.csv');wr(bind(worstV),'maximum_V_share_examples.csv')
 # Independent source-headline checks,with distinct cohort/metric finite masks.
 for(st in c('architecture','stress')){
  p<-file.path(P,'Scientific Reports Manuscript/github_results/summary_results',if(st=='architecture')'architecture_M1_M4_pooled.csv' else 'stress_pooled.csv');z<-read.csv(p);z<-z[z$method=='OTRIMLE_fixedK',];stopifnot(nrow(z)==1)
  co<-if(st=='architecture')'Original architecture M1-M4' else 'Original stress';sm<-bind(cohort_rates)
  for(m in c('MajorARI','MinorityClusterF1','MinorityDetectionF1')){
   actual<-sm$mean[sm$cohort==co&sm$metric==m];stopifnot(length(actual)==1,abs(actual-z[[m]])<1e-11)
   headline[[length(headline)+1]]<<-data.frame(cohort=co,metric=m,n_fits=z$N,n_finite=sm$n_finite[sm$cohort==co&sm$metric==m],reproduced_mean=actual,submitted_headline=z[[m]],absolute_difference=abs(actual-z[[m]]),tolerance=1e-11,source=p)
  }
 }
 wr(bind(headline),'headline_checks.csv','evidence')
 wr(data.frame(check=c('235 immutable scientific sources','54200 original OTRIMLE records and all66 exported fields','400 independent-context D1-E OTRIMLE rows','Independent integer-flow versus analytic bounds','All prior per-fit aggregate counts and endpoints','Nested promotion inequalities and endpoint flow witnesses','Metric-specific finite masks and original inclusion rules'),status='PASS',detail=c('SHA256 checked before analysis','RDS and authoritative gzip exports agree','No original/revision pooling;all stored macroE rates0',paste(witness_n,'records;previous bound function never called'),'Exact integer/count equality for every keyed saved fit','Nonnegative integer flow witnesses at both V endpoints;relaxation not reconstructed labels','Negative-control minority scores undefined;code1 retained')),'analysis_checks.csv','evidence')
 writeLines(capture.output(sessionInfo()),file.path(F,'evidence/R_session_info.txt'))
 cat('INDEPENDENT_FINAL_BOUNDS_OK',nrow(all),'records;NO_FITS\n');print(bind(cohort_summary)[,c('cohort','target','n_fits','n_identified','n_nonunique','sum_lower_count','sum_upper_count','mean_lower_pct','mean_upper_pct','max_individual_upper_pct')]);print(bind(collapse))
}
main()
